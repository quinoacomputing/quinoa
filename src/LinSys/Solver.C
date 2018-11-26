// *****************************************************************************
/*!
  \file      src/LinSys/Solver.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ linear system merger nodegroup to solve a linear system
  \see       src/LinSys/Solver.h
*/
// *****************************************************************************

#include <numeric>

#include "Exception.h"
#include "ContainerUtil.h"
#include "VectorReducer.h"
#include "HashMapReducer.h"

#include "Solver.h"

using tk::Solver;

Solver::Solver( const SolverCallback& cb, std::size_t n ) :
  m_cb( cb ),
  m_ncomp( n ),
  m_nchare( 0 ),
  m_mynchare( 0 ),
  m_nbounds( 0 ),
  m_nchbc( 0 ),
  m_lower( std::numeric_limits< std::size_t >::max() ),
  m_upper( 0 ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_dt( 0.0 ),
  m_initial( true ),
  m_worker(),
  m_rowimport(),
  m_solimport(),
  m_lhsimport(),
  m_rhsimport(),
  m_lowrhsimport(),
  m_lowlhsimport(),
  m_row(),
  m_sol(),
  m_lhs(),
  m_rhs(),
  m_lowrhs(),
  m_lowlhs(),
  m_x(),
  m_A(),
  m_b(),
  m_solver(),
  m_hypreRows(),
  m_hypreNcols(),
  m_hypreCols(),
  m_hypreMat(),
  m_hypreRhs(),
  m_hypreSol(),
  m_lid(),
  m_div(),
  m_node(),
  m_bc(),
  m_bca()
// *****************************************************************************
//  Constructor
//! \param[in] cb Charm++ callbacks for Solver
//! \param[in] n Total number of scalar components in the linear system
// *****************************************************************************
{
  // Activate SDAG waits
  thisProxy[ CkMyNode() ].wait4com();
  thisProxy[ CkMyNode() ].wait4lhsbc();
  thisProxy[ CkMyNode() ].wait4rhsbc();
  thisProxy[ CkMyNode() ].wait4hypresol();
  thisProxy[ CkMyNode() ].wait4hyprelhs();
  thisProxy[ CkMyNode() ].wait4hyprerhs();
  thisProxy[ CkMyNode() ].wait4asm();
  thisProxy[ CkMyNode() ].wait4low();
}

void
Solver::nchare( int n )
// *****************************************************************************
//  Set number of worker chares expected to contribute on this compute node
//! \param[in] n Total number of chares (work units) across all compute nodes
//! \details Besides the number of workers contribute to this compute node, we
//!    also store the total number of chares across the whole problem.
// *****************************************************************************
{
  auto chunksize = n / CkNumNodes();
  auto mynchare = chunksize;
  if (CkMyNode() == CkNumNodes()-1) mynchare += n % CkNumNodes();
  m_mynchare = static_cast< std::size_t >( mynchare );
  m_nchare = static_cast< std::size_t >( n );

  contribute( m_cb.get< tag::part >() );
}

void
Solver::chbounds( std::size_t lower, std::size_t upper )
// *****************************************************************************
//  Receive lower and upper global node IDs from chares
//! \param[in] lower Lower index of the global rows of sending chare
//! \param[in] upper Upper index of the global rows of sending chare
// *****************************************************************************
{
  Assert( lower < upper, "Lower bound must be lower than the upper bound: "
          "(" + std::to_string(lower) + "..." +  std::to_string(upper) +
          ") sent by chare" );

  // Compute extents of bounds on this compute node
  m_lower = std::min( m_lower, lower );
  m_upper = std::max( m_upper, upper );

  // When this compute node has heard from all chares on this compute node,
  // aggregate compute node bounds
  if (++m_nbounds == m_mynchare)
    thisProxy.nodebounds( CkMyNode(), m_lower, m_upper );
}

void
Solver::nodebounds( int n, std::size_t lower, std::size_t upper )
// *****************************************************************************
//  Receive lower and upper bounds across all compute nodes
//! \param[in] n Compute node whose bounds being received
//! \param[in] lower Lower index of the global rows of sending compute node
//! \param[in] upper Upper index of the global rows of sending compute node
// *****************************************************************************
{
  Assert( lower < upper, "Lower bound must be lower than the upper bound: "
          "(" + std::to_string(lower) + "..." +  std::to_string(upper) +
          ") sent by compute node " + std::to_string(n) );

  // Store our bounds
  if (n == CkMyNode()) {
    m_lower = lower;
    m_upper = upper;
  }

  // Store inverse of compute-node-division map stored on all compute nodes
  m_div[ {lower,upper} ] = n;

  // If we have all compute nodes' bounds, signal the runtime system to continue
  if (m_div.size() == static_cast<std::size_t>(CkNumNodes())) {

    // Create my compute node's lhs matrix distributed across all compute nodes
    m_A.create( m_lower*m_ncomp, m_upper*m_ncomp );
    // Create my compute node's rhs and unknown vectors distributed across all
    // compute nodes
    m_b.create( m_lower*m_ncomp, m_upper*m_ncomp );
    m_x.create( m_lower*m_ncomp, m_upper*m_ncomp );
    // Create linear solver
    m_solver.create();
    contribute( m_cb.get< tag::bounds >() );
  }
}

void
Solver::next()
// *****************************************************************************
//  Prepare for next step
//! \details Re-enable SDAG waits for rebuilding the right-hand side vector only
// *****************************************************************************
{
  m_initial = false;

  thisProxy[ CkMyNode() ].wait4rhsbc();
  thisProxy[ CkMyNode() ].wait4hyprerhs();
  thisProxy[ CkMyNode() ].wait4asm();
  thisProxy[ CkMyNode() ].wait4low();

  m_rhsimport.clear();
  m_lowrhsimport.clear();
  m_lowrhs.clear();
  m_hypreRhs.clear();
  m_rhs.clear();

  m_bc.clear();

  lowlhs_complete();
  hyprerow_complete();
  asmsol_complete();
  asmlhs_complete();

  // Continue with next time step: for all workers call .dt()
  if (CkMyNode() == 0)
    for (const auto& c : m_worker)
       c.second.get< tag::dt >().send();
}

void
Solver::charecom( int fromch, const MatCGCallback& cb )
// *****************************************************************************
//  Chares contribute their ids and callbacks
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] cb Callbacks to member functions
// *****************************************************************************
{
  // Store ids and callbacks to worker chare entry methods
  m_worker[ fromch ] = cb;

  if (m_worker.size() == m_nchare) com_complete();
}

void
Solver::charerow( int fromch, const std::vector< std::size_t >& row )
// *****************************************************************************
//  Chares contribute their global row ids for establishing communications
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] row Global mesh point (row) indices contributed
// *****************************************************************************
{
  // Store rows owned and pack those to be exported, also build import map
  // used to test for completion
  std::map< int, std::set< std::size_t > > exp;
  for (auto gid : row) {
    if (gid >= m_lower && gid < m_upper) {  // if own
      m_rowimport[ fromch ].push_back( gid );
      m_row.insert( gid );
    } else exp[ node(gid) ].insert( gid );
  }

  // Export non-owned parts to fellow branches that own them
  for (const auto& p : exp) {
    auto tonode = static_cast< int >( p.first );
    thisProxy[ tonode ].addrow( fromch, p.second );
  }

  if (m_row.size() == m_upper-m_lower) row_complete();
}

void
Solver::addrow( int fromch, const std::set< std::size_t >& row )
// *****************************************************************************
//  Receive global row ids from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] row Global mesh point (row) indices received
// *****************************************************************************
{
  for (auto r : row) {
    m_rowimport[ fromch ].push_back( r );
    m_row.insert( r );
  }

  if (m_row.size() == m_upper-m_lower) row_complete();
}

void
Solver::charesol( int fromch,
                  const std::vector< std::size_t >& gid,
                  const Fields& solution )
// *****************************************************************************
//  Chares contribute their solution nonzero values
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] solution Portion of the unknown/solution vector contributed
// *****************************************************************************
{
  Assert( gid.size() == solution.nunk(),
          "Size of solution and row ID vectors must equal" );

  // Store solution vector nonzero values owned and pack those to be
  // exported, also build import map used to test for completion
  std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;

  for (std::size_t i=0; i<gid.size(); ++i)
    if (gid[i] >= m_lower && gid[i] < m_upper) {    // if own
      m_solimport[ fromch ].push_back( gid[i] );
      m_sol[ gid[i] ] = solution[i];
    } else {
      exp[ node(gid[i]) ][ gid[i] ] = solution[i];
    }

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tonode = static_cast< int >( p.first );
    thisProxy[ tonode ].addsol( fromch, p.second );
  }

  if (solcomplete()) hypresol();
}

void
Solver::addsol( int fromch,
                const std::map< std::size_t,
                                std::vector< tk::real > >& solution )
// *****************************************************************************
//  Receive solution vector nonzeros from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] solution Portion of the unknown/solution vector contributed,
//!   containing global row indices and values for all components
// *****************************************************************************
{
  for (const auto& r : solution) {
    m_solimport[ fromch ].push_back( r.first );
    m_sol[ r.first ] = r.second;
  }

  if (solcomplete()) hypresol();
}

void
Solver::charelhs( int fromch,
                  const std::vector< std::size_t >& gid,
                  const std::pair< std::vector< std::size_t >,
                                   std::vector< std::size_t > >& psup,
                  const tk::Fields& lhsd,
                  const tk::Fields& lhso )
// *****************************************************************************
//  Chares contribute their matrix nonzero values
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] gid Global row indices of the matrix chunk contributed
//! \param[in] psup Points surrounding points using local indices. See also
//!   tk::genPsup().
//! \param[in] lhsd Portion of the left-hand side matrix contributed,
//!   containing non-zero values (for all scalar components of the equations
//!   solved) as a sparse matrix diagonal
//! \param[in] lhso Portion of the left-hand side matrix contributed,
//!   containing non-zero values (for all scalar components of the equations
//!   solved) as a sparse matrix off-diagonal entries in compressed row
//!   storage format
// *****************************************************************************
{
  Assert( psup.second.size()-1 == gid.size(),
          "Number of mesh points and number of global IDs unequal" );
  Assert( psup.second.size()-1 == lhsd.nunk(),
          "Number of mesh points and number of diagonals unequal" );
  Assert( psup.first.size() == lhso.nunk(),
          "Number of off-diagonals and their number of indices unequal" );

  // Store matrix nonzero values owned and pack those to be exported, also
  // build import map used to test for completion
  std::map< int, std::map< std::size_t,
                           std::map< std::size_t,
                                     std::vector< tk::real > > > > exp;

  for (std::size_t i=0; i<gid.size(); ++i)
    if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
      m_lhsimport[ fromch ].push_back( gid[i] );
      auto& row = m_lhs[ gid[i] ];
      row[ gid[i] ] += lhsd[i];
      for (auto j=psup.second[i]+1; j<=psup.second[i+1]; ++j)
        row[ gid[ psup.first[j] ] ] += lhso[j];
    } else {
      auto& row = exp[ node(gid[i]) ][ gid[i] ];
      row[ gid[i] ] = lhsd[i];
      for (auto j=psup.second[i]+1; j<=psup.second[i+1]; ++j)
        row[ gid[ psup.first[j] ] ] = lhso[j];
    }

  // Export non-owned matrix rows values to fellow branches that own them
  for (const auto& p : exp)
    thisProxy[ p.first ].addlhs( fromch, p.second );

  if (lhscomplete()) lhs_complete();
}

void
Solver::addlhs( int fromch,
                const std::map< std::size_t,
                                std::map< std::size_t,
                                          std::vector< tk::real > > >& l )
// *****************************************************************************
//  Receive matrix nonzeros from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] l Portion of the left-hand side matrix contributed,
//!   containing global row and column indices and non-zero values
// *****************************************************************************
{
  for (const auto& r : l) {
    m_lhsimport[ fromch ].push_back( r.first );
    auto& row = m_lhs[ r.first ];
    for (const auto& c : r.second) row[ c.first ] += c.second;
  }

  if (lhscomplete()) lhs_complete();
}

void
Solver::charerhs( int fromch,
                  const std::vector< std::size_t >& gid,
                  const Fields& r )
// *****************************************************************************
//  Chares contribute their rhs nonzero values
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] r Portion of the right-hand side vector contributed
// *****************************************************************************
{
  Assert( gid.size() == r.nunk(),
          "Size of right-hand side and row ID vectors must equal" );

  // Store+add vector of nonzero values owned and pack those to be exported
  std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;

  for (std::size_t i=0; i<gid.size(); ++i)
    if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
      m_rhsimport[ fromch ].push_back( gid[i] );
      m_rhs[ gid[i] ] += r[i];
    } else
      exp[ node(gid[i]) ][ gid[i] ] = r[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tonode = static_cast< int >( p.first );
    thisProxy[ tonode ].addrhs( fromch, p.second );
  }

  if (rhscomplete()) rhs_complete();
}

void
Solver::addrhs( int fromch,
                const std::map< std::size_t, std::vector< tk::real > >& r )
// *****************************************************************************
//  Receive+add right-hand side vector nonzeros from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] r Portion of the right-hand side vector contributed, containing
//!   global row indices and values
// *****************************************************************************
{
  // Store rhs contributions
  for (const auto& l : r) {
    m_rhsimport[ fromch ].push_back( l.first );
    m_rhs[ l.first ] += l.second;
  }

  if (rhscomplete()) rhs_complete();
}

void
Solver::charelowrhs( int fromch,
                     const std::vector< std::size_t >& gid,
                     const Fields& lowrhs )
// *****************************************************************************
//  Chares contribute to the rhs of the low order linear system
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] lowrhs Portion of the low order rhs vector contributed
// *****************************************************************************
{
  using tk::operator+=;

  Assert( gid.size() == lowrhs.nunk(),
          "Size of mass diffusion rhs and row ID vectors must equal" );

  // Store+add vector of nonzero values owned and pack those to be exported
  std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;

  for (std::size_t i=0; i<gid.size(); ++i)
    if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
      m_lowrhsimport[ fromch ].push_back( gid[i] );
      m_lowrhs[ gid[i] ] += lowrhs[i];
    } else
      exp[ node(gid[i]) ][ gid[i] ] = lowrhs[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tonode = static_cast< int >( p.first );
    thisProxy[ tonode ].addlowrhs( fromch, p.second );
  }

  if (lowrhscomplete()) lowrhs_complete();
}

void
Solver::addlowrhs( int fromch,
                   const std::map< std::size_t,
                                   std::vector< tk::real > >& lowrhs )
// *****************************************************************************
//  Receive+add low order rhs vector nonzeros from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] lowrhs Portion of the low order rhs vector contributed,
//!   containing global row indices and values
// *****************************************************************************
{
  using tk::operator+=;

  for (const auto& r : lowrhs) {
    m_lowrhsimport[ fromch ].push_back( r.first );
    m_lowrhs[ r.first ] += r.second;
  }

  if (lowrhscomplete()) lowrhs_complete();
}

void
Solver::charelowlhs( int fromch,
                     const std::vector< std::size_t >& gid,
                     const Fields& lowlhs )
// *****************************************************************************
//  Chares contribute to the lhs of the low order linear system
//! \param[in] fromch Charm chare array the contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] lowlhs Portion of the low order lhs vector contributed
// *****************************************************************************
{
  using tk::operator+=;

  Assert( gid.size() == lowlhs.nunk(),
          "Size of mass diffusion lhs and row ID vectors must equal" );

  // Store+add vector of nonzero values owned and pack those to be exported
  std::map< int, std::map< std::size_t, std::vector< tk::real > > > exp;

  for (std::size_t i=0; i<gid.size(); ++i)
    if (gid[i] >= m_lower && gid[i] < m_upper) {  // if own
      m_lowlhsimport[ fromch ].push_back( gid[i] );
      m_lowlhs[ gid[i] ] += lowlhs[i];
    } else
      exp[ node(gid[i]) ][ gid[i] ] = lowlhs[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tonode = static_cast< int >( p.first );
    thisProxy[ tonode ].addlowlhs( fromch, p.second );
  }

  if (lowlhscomplete()) lowlhs_complete();
}

void
Solver::addlowlhs( int fromch,
                   const std::map< std::size_t,
                                   std::vector< tk::real > >& lowlhs )
// *****************************************************************************
// Receive and add lhs vector to the low order system from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] lowlhs Portion of the lhs vector contributed to the low order
//!   linear system, containing global row indices and values
// *****************************************************************************
{
  using tk::operator+=;

  for (const auto& r : lowlhs) {
    m_lowlhsimport[ fromch ].push_back( r.first );
    m_lowlhs[ r.first ] += r.second;
  }

  if (lowlhscomplete()) lowlhs_complete();
}

void
Solver::comfinal()
// *****************************************************************************
//  All communications have been establised among compute nodes
//! \details At this point all solver objects on all compute nodes must have
//!   received their global row ids which means that the communications have
//!   been established among all compute nodes and this the communications
//!   (maps) are final on all compute nodes.
// *****************************************************************************
{
  Assert( m_row.size() == m_upper-m_lower,
          "Row ids are incomplete on node " + std::to_string(CkMyNode()) + ": "
          "number of rows received: " + std::to_string(m_row.size()) + " vs. "
          "number of rows supposed to have been received: " +
          std::to_string(m_upper-m_lower) );

  // now that the global row ids are complete, build Hypre data from it
  hyprerow();
}

void
Solver::charebc( const std::unordered_map< std::size_t,
                   std::vector< std::pair< bool, tk::real > > >& bc )
// *****************************************************************************
//  Chares contribute their global row ids and associated Dirichlet boundary
//  condition values at which they set BCs
//! \param[in] bc Vector of pairs of bool and BC value (BC vector)
//!   associated to global node IDs at which the boundary condition is set.
//!   Here the bool indicates whether the BC value is set at the given node
//!   by the user. The size of the vectors is the number of PDEs integrated
//!   times the number of scalar components in all PDEs.
// *****************************************************************************
{
  // Associate BC vectors to mesh nodes owned
  for (const auto& n : bc) {
    Assert( n.second.size() == m_ncomp, "The total number of scalar "
            "components does not equal that of set in the BC vector." );
    m_bc[ n.first ] = n.second;
  }

  if (++m_nchbc == m_nchare) {
    m_nchbc = 0;
    bc_complete();
    if (m_initial) bc_complete();
  }
}

int
Solver::node( std::size_t gid )
// *****************************************************************************
//  Return compute node id for global mesh row id
//! \param[in] gid Global mesh point (matrix or vector row) id
//! \details First we attempt to find the point index in the cache. If that
//!   fails, we resort to a linear search across the division map. Once the
//!   compute node is found, we store it in the cache, so next time the search
//!   is quicker. This procedure must find the compute node id for the global
//!   node/row id.
//! \return Compute node id that owns the global row id
// *****************************************************************************
{
  Assert( m_div.size() == static_cast<std::size_t>(CkNumNodes()),
          "Compute node bounds incomplete on node " +
          std::to_string(CkMyNode()) );

  int n = -1;
  auto it = m_node.find( gid );
  if (it != end(m_node))
    n = it->second;
  else
    for (const auto& d : m_div)
      if (gid >= d.first.first && gid < d.first.second)
        n = m_node[ gid ] = d.second;

  Assert( n >= 0, "Compute node not found for node id " + std::to_string(gid) );
  Assert( n < CkNumNodes(), "Assigning to nonexistent compute node" );

  return n;
}

void
Solver::hyprerow()
// *****************************************************************************
//  Build Hypre data for our portion of the global row ids
//! \note Hypre only likes one-based indexing. Zero-based row indexing fails
//!   to update the vector with HYPRE_IJVectorGetValues().
// *****************************************************************************
{
  if (m_hypreRows.empty()) {

    for (auto r : m_row) {
      std::vector< int > h( m_ncomp );
      std::iota( begin(h), end(h), r*m_ncomp+1 );
      m_hypreRows.insert( end(m_hypreRows), begin(h), end(h) );
    }

    hyprerow_complete();
    if (m_initial) {
      hyprerow_complete();
      hyprerow_complete();
    }

  }
}

void
Solver::lhsbc()
// *****************************************************************************
//  Set boundary conditions on the left-hand side matrix
//! \details Setting boundary conditions on the left-hand side matrix is
//!   done by zeroing all nonzero entries in rows where BCs are prescribed
//!   and putting in 1.0 for the diagonal.
// *****************************************************************************
{
  Assert( lhscomplete(),
          "Nonzero values of distributed matrix on compute node " +
          std::to_string( CkMyNode() ) + " is incomplete: cannot set BCs" );

  // Set Dirichlet BCs on the lhs matrix. Loop through all BCs and if a BC
  // is prescribed on a row we own, find that row (r) and in that row the
  // position of the diagonal entry (diag) in the LHS matrix. Then for all
  // scalar components the system of system of PDEs we intagrate, query the
  // entry in the BC data structure to see if we need to set BC for the
  // given component and set the off-diagonals to zero while put 1.0 into
  // the diagonal.

//std::cout << "l: ";
// std::map< std::size_t, std::vector< std::pair< bool, tk::real > > >
//   s( begin(m_bc), end(m_bc) );
// for (const auto& n : s) std::cout << n.first << ' '; std::cout << '\n';

  for (const auto& n : m_bc) {
    if (n.first >= m_lower && n.first < m_upper) {
      auto& r = tk::ref_find( m_lhs, n.first );
      auto& diag = tk::ref_find( r, n.first );
      for (std::size_t i=0; i<m_ncomp; ++i)
        if (n.second[i].first) {
          for (auto& c : r) c.second[i] = 0.0;  // zero columns in BC row
          diag[i] = 1.0;    // put 1.0 in diagonal of BC row
        }
    }
  }

//   for (const auto& n : m_bc) {
//     auto r = m_lhs.find( n.first );
//     if (r != end(m_lhs)) {
//       m_bca[ n.first ] = r->second;
//       for (const auto& c : r->second)
//         if (c.first != n.first ) {
//           auto a = m_lhs.find( c.first );
//           if (a != end(m_lhs)) {
//             auto& b = tk::ref_find( a->second, n.first );
//             for (std::size_t i=0; i<m_ncomp; ++i)
//               if (n.second[i].first)
//                 b[i] = 0.0;
//           }
//         }
//     }
//   }
// 
//   // test BCs in matrix (only in serial)
//   auto eps = std::numeric_limits< tk::real >::epsilon();
//   for (const auto& n : m_bc) {
//     auto r = m_lhs.find( n.first );
//     if (r != end(m_lhs)) {
//       std::size_t ver_row = 0, ver_col = 0;
//       for (const auto& c : r->second)
//         if (c.first != n.first) {
//           // test row
//           for (std::size_t i=0; i<m_ncomp; ++i)
//             if (std::abs(c.second[i]) > eps)      // test zero
//               std::cout << 'r';
//             else                                  // count nonzeros
//               ++ver_row;
//           // test column
//           auto b = m_lhs.find( c.first );
//           if (b !=end(m_lhs)) {
//             auto d = b->second.find( n.first );
//             if (d != end(b->second)) {
//               for (std::size_t i=0; i<m_ncomp; ++i)
//                 if (std::abs(d->second[i]) > eps)  // test zero
//                   std::cout << 'c';
//                 else
//                   ++ver_col;                       // count nonzeros
//             }
//           }
//         }
//       // test number of nonzeros in BC row
//       if (ver_row != (r->second.size()-1)*m_ncomp) std::cout << 'R';
//       // test number of nonzeros in BC column
//       if (ver_col != (r->second.size()-1)*m_ncomp) std::cout << 'C';
//     }
//   }
// 
//   // test symmetry of matrix (only in serial)
//   for (const auto& r : m_lhs)
//     for (const auto& c : r.second) {
//       const auto& v = tk::cref_find( m_lhs, c.first );
//       const auto& w = tk::cref_find( v, r.first );
//       for (std::size_t i=0; i<m_ncomp; ++i)
//         if (std::abs(c.second[i]-w[i]) > 1.0e-14)
//           std::cout << 's';
//     }

  hyprelhs();
}

void
Solver::rhsbc()
// *****************************************************************************
//  Set Dirichlet boundary conditions on the right-hand side vector
//! \details Since we solve for the solution increment, this amounts to
//!    enforcing zero rhs (no solution increment) at BC nodes
// *****************************************************************************
{
  Assert( rhscomplete(), "Values of distributed right-hand-side vector on "
          "compute node " + std::to_string( CkMyNode() ) + " is incomplete: "
          "cannot set BCs" );

//   for (const auto& n : m_bc) {
//      auto r = m_bca.find( n.first );
//      if (r != end(m_bca))
//        for (const auto& c : r->second)
//          if (c.first != n.first ) {
//            auto a = m_rhs.find( c.first );
//            if (a != end(m_rhs))
//              for (std::size_t i=0; i<m_ncomp; ++i)
//                if (n.second[i].first)
//                  a->second[i] -= n.second[i].second * c.second[i];
//          }
//    }

  for (const auto& n : m_bc) {
    if (n.first >= m_lower && n.first < m_upper) {
      auto& r = tk::ref_find( m_rhs, n.first );
      for (std::size_t i=0; i<m_ncomp; ++i)
        if (n.second[i].first)
          r[i] = n.second[i].second;
    }
  }

  hyprerhs();

  rhsbc_complete();
}

void
Solver::hypresol()
// *****************************************************************************
//  Build Hypre data for our portion of the solution vector
// *****************************************************************************
{
  Assert( solcomplete(), "Values of distributed solution vector on compute "
          "node " + std::to_string( CkMyNode() ) + " is incomplete" );

  std::size_t i = 0;
  for (const auto& r : m_sol) {
    m_lid[ r.first ] = i++;
    m_hypreSol.insert( end(m_hypreSol), begin(r.second), end(r.second) );
  }

  hypresol_complete();
}

void
Solver::hyprelhs()
// *****************************************************************************
//  Build Hypre data for our portion of the matrix
//! \note Hypre only likes one-based indexing. Zero-based row indexing fails
//!   to update the vector with HYPRE_IJVectorGetValues().
// *****************************************************************************
{
  Assert( lhscomplete(), "Nonzero values of distributed matrix on compute "
          "node " + std::to_string( CkMyNode() ) + " is incomplete: cannot "
          "convert" );

  for (const auto& r : m_lhs)
    for (std::size_t i=0; i<m_ncomp; ++i) {
      m_hypreNcols.push_back( static_cast< int >( r.second.size() ) );
      for (const auto& c : r.second) {
        m_hypreCols.push_back( static_cast< int >( c.first*m_ncomp+i+1 ) );
        m_hypreMat.push_back( c.second[i] );
      }
    }

  hyprelhs_complete();
}

void
Solver::hyprerhs()
// *****************************************************************************
//  Build Hypre data for our portion of the right-hand side vector
// *****************************************************************************
{
  Assert( rhscomplete(), "Values of distributed right-hand-side vector on "
          "compute node " + std::to_string( CkMyNode() ) + " is incomplete: "
          "cannot convert" );

  for (const auto& r : m_rhs)
    m_hypreRhs.insert( end(m_hypreRhs), begin(r.second), end(r.second) );

  hyprerhs_complete();
}

void
Solver::sol()
// *****************************************************************************
//  Set our portion of values of the distributed solution vector
// *****************************************************************************
{
  Assert( m_hypreSol.size() == m_hypreRows.size(), "Solution vector values "
          "incomplete on compute node " + std::to_string(CkMyNode()) );

  // Set our portion of the vector values
  m_x.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
           m_hypreRows.data(),
           m_hypreSol.data() );

  m_x.assemble();
  asmsol_complete();
}

void
Solver::lhs()
// *****************************************************************************
//  Set our portion of values of the distributed matrix
// *****************************************************************************
{
  Assert( m_hypreMat.size() == m_hypreCols.size(), "Matrix values incomplete "
          "on compute node " + std::to_string(CkMyNode()) );

  // Set our portion of the matrix values
  m_A.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
           m_hypreNcols.data(),
           m_hypreRows.data(),
           m_hypreCols.data(),
           m_hypreMat.data() );

  m_A.assemble();
  asmlhs_complete();
}

void
Solver::rhs()
// *****************************************************************************
//  Set our portion of values of the distributed right-hand side vector
// *****************************************************************************
{
  Assert( m_hypreRhs.size() == m_hypreRows.size(), "RHS vector values "
          "incomplete on compute node " + std::to_string(CkMyNode()) );

  // Set our portion of the vector values
  m_b.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
           m_hypreRows.data(),
           m_hypreRhs.data() );

  m_b.assemble();
  asmrhs_complete();
}

std::pair< int, std::unique_ptr<char[]> >
Solver::serializeSol( const std::vector< std::size_t >& gid,
                      const std::vector< tk::real >& u ) const
// *****************************************************************************
//  Serialize solution vector into a Charm++ message, ready for a CkCallback
//! \param[in] gid Global row indices of the vector updated
//! \param[in] u Portion of the unknown/solution vector update
// *****************************************************************************
{
  // Prepare for serializing global ids and solution vector to a raw binary
  // stream, compute size
  PUP::sizer sizer;
  sizer | const_cast< std::vector< std::size_t >& >( gid );
  sizer | const_cast< std::vector< tk::real >& >( u );

  // Create raw character stream to store the serialized data
  std::unique_ptr<char[]> flatData = tk::make_unique<char[]>( sizer.size() );

  // Serialize global ids and solution
  PUP::toMem packer( flatData.get() );
  packer | const_cast< std::vector< std::size_t >& >( gid );
  packer | const_cast< std::vector< tk::real >& >( u );

  // Return serialized data (size in number of characters and pointer to data)
  return { sizer.size(), std::move(flatData) };
}

void
Solver::updateSol()
// *****************************************************************************
//  Update solution vector in workers on this compute node
// *****************************************************************************
{
  // Get solution vector values for this compute node
  m_x.get( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
           m_hypreRows.data(),
           m_hypreSol.data() );

  // Group solution vector by workers and send each the parts back to workers
  // that own them
  for (const auto& w : m_solimport) {
    std::vector< std::size_t > gid;
    std::vector< tk::real > solution;

    for (auto r : w.second) {
      const auto it = m_sol.find( r );
      if (it != end(m_sol)) {
        gid.push_back( it->first );
        auto i = tk::cref_find( m_lid, it->first );
        using diff_type = typename decltype(m_hypreSol)::difference_type;
        auto b = static_cast< diff_type >( i*m_ncomp );
        auto e = static_cast< diff_type >( (i+1)*m_ncomp );
        solution.insert( end(solution),
                         std::next( begin(m_hypreSol), b ),
                         std::next( begin(m_hypreSol), e ) );
      } else
        Throw( "Can't find global row id " + std::to_string(r) +
               " to export in solution vector" );
    }

    // Update worker with high order solution
    auto stream = serializeSol( gid, solution );
    tk::cref_find( m_worker, w.first ).get< tag::high >().
      send( stream.first, stream.second.get() );
  }
}

void
Solver::solve()
// *****************************************************************************
//  Solve high order linear system
// *****************************************************************************
{
  m_solver.solve( m_A, m_b, m_x );
  updateSol();
}

void
Solver::updateLowSol()
// *****************************************************************************
//  Update low order solution vector in workers on this compute node
// *****************************************************************************
{
  for (const auto& w : m_solimport) {
    std::vector< std::size_t > gid;
    std::vector< tk::real > solution;

    for (auto r : w.second) {
      const auto it = m_lowrhs.find( r );
      if (it != end(m_lowrhs)) {
        gid.push_back( it->first );
        solution.insert(end(solution), begin(it->second), end(it->second));
      } else
        Throw( "Can't find global row id " + std::to_string(r) +
               " to export in low order solution vector" );
    }

    // Update worker with low order solution
    auto stream = serializeSol( gid, solution );
    tk::cref_find( m_worker, w.first ).get< tag::low >().
      send( stream.first, stream.second.get() );
  }
}

void
Solver::lowsolve()
// *****************************************************************************
//  Solve low order linear system
// *****************************************************************************
{
  // Set boundary conditions on the low order system
  for (const auto& n : m_bc)
    if (n.first >= m_lower && n.first < m_upper) {
      // lhs
      auto& l = tk::ref_find( m_lowlhs, n.first );
      for (std::size_t i=0; i<m_ncomp; ++i)
        if (n.second[i].first) l[i] = 1.0;
      // rhs (set to zero instead of the solution increment at Dirichlet
      // BCs, because for the low order solution the right hand side is the sum
      // of the high order right hand side and mass diffusion, so the low order
      // system is L = R + D, where L is the lumped mass matrix, R is the high
      // order RHS, and D is mass diffusion, and R already has the Dirichlet BC
      // set)
      auto& r = tk::ref_find( m_lowrhs, n.first );
      for (std::size_t i=0; i<m_ncomp; ++i)
        if (n.second[i].first) r[i] = 0.0;
    }

  // Solve low order system
  Assert( rhscomplete(),
          "Values of distributed right-hand-side vector on compute node " +
          std::to_string( CkMyNode() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( lowrhscomplete(),
          "Values of distributed mass diffusion rhs vector on compute node " +
          std::to_string( CkMyNode() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( lowlhscomplete(),
          "Values of distributed lumped mass lhs vector on compute node " +
          std::to_string( CkMyNode() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( tk::keyEqual( m_rhs, m_lowrhs ), "Row IDs of rhs and mass "
          "diffusion rhs vector unequal on compute node " +
          std::to_string( CkMyNode() )  + ": cannot solve low order system" );
  Assert( tk::keyEqual( m_rhs, m_lowlhs ), "Row IDs of rhs and lumped mass "
          "lhs vector unequal on compute node " + std::to_string( CkMyNode() ) +
          ": cannot solve low order system" );

  auto ir = m_rhs.cbegin();
  auto id = m_lowrhs.begin();
  auto im = m_lowlhs.cbegin();

  while (ir != m_rhs.cend()) {
    const auto& r = ir->second;
    const auto& m = im->second;
    auto& d = id->second;
    Assert( r.size()==m_ncomp && m.size()==m_ncomp && d.size()==m_ncomp,
            "Wrong number of components in solving the low order system" );
    for (std::size_t i=0; i<m_ncomp; ++i) d[i] = (r[i]+d[i])/m[i];
    ++ir; ++id; ++im;
  }

  //lowsolve_complete();
  updateLowSol();
}

#include "NoWarning/solver.def.h"
