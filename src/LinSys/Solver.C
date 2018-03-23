// *****************************************************************************
/*!
  \file      src/LinSys/Solver.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Charm++ chare linear system merger group to solve a linear system
  \details   Charm++ chare linear system merger group used to collect and
    assemble the left hand side matrix (lhs), the right hand side (rhs) vector,
    and the solution (unknown) vector from individual worker
    chares. Beside collection and assembly, the system is also solved. The
    solution is outsourced to hypre, an MPI-only library. Once the solution is
    available, the individual worker chares are updated with the new solution.

    The implementation uses the Charm++ runtime system and is fully
    asynchronous, overlapping computation and communication. The algorithm
    utilizes the structured dagger (SDAG) Charm++ functionality. The high-level
    overview of the algorithm structure and how it interfaces with Charm++ is
    discussed in the Charm++ interface file src/LinSys/solver.ci. See also
    src/LinSys/Solver.h for a discussion of the asynchronous call logic.
*/
// *****************************************************************************

#include <numeric>

#include "Exception.h"
#include "ContainerUtil.h"
#include "VectorReducer.h"
#include "HashMapReducer.h"

#include "Solver.h"

namespace tk {

static CkReduction::reducerType BCMapMerger;

}

tk::SolverShadow::SolverShadow()
// *****************************************************************************
//  Constructor
//! \details Solver shadow class constructor used to fire off a reduction
//!   different from Solver to avoid the runtime error "mis-matched client
//!   callbacks in reduction messages". Why the constructor definition is not
//!   defined in the class definition? To avoid the compiler warning: "warning:
//!   instantiation of function 'CBaseT1<Group,
//!   tk::CProxy_SolverShadow>::virtual_pup' required here, but no definition is
//!   available [-Wundefined-func-template]".
// *****************************************************************************
{
}

using tk::Solver;

Solver::Solver( CProxy_SolverShadow sh,
                const std::vector< CkCallback >& cb,
                std::size_t n,
                bool /*feedback*/ ) :
  m_shadow( sh ),
  m_cb( cb[0], cb[1], cb[2] ),
  m_ncomp( n ),
  m_nchare( 0 ),
  m_ncomm( 0 ),
  m_nperow( 0 ),
  m_nchbc( 0 ),
  m_lower( 0 ),
  m_upper( 0 ),
  m_it( 0 ),
  m_t( 0.0 ),
  m_dt( 0.0 ),
  //m_feedback( feedback ),
  m_myworker(),
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
  m_pe(),
  m_bc(),
  m_bca()
// *****************************************************************************
//  Constructor
//! \param[in] cb Charm++ callbacks
//! \param[in] s Mesh node IDs mapped to side set ids
//! \param[in] n Total number of scalar components in the linear system
// *****************************************************************************
{
  // Activate SDAG waits
  wait4nchare();
  wait4lhsbc();
  wait4rhsbc();
  wait4hypresol();
  wait4hyprelhs();
  wait4hyprerhs();
  wait4asm();
  wait4low();
}

void
Solver::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  BCMapMerger = CkReduction::addReducer(
                  tk::mergeHashMap< std::size_t,
                    std::vector< std::pair< bool, tk::real > > > );
}

void
Solver::bounds( int p, std::size_t lower, std::size_t upper )
// *****************************************************************************
//  Receive lower and upper global node IDs all PEs will operate on
//! \param[in] p PE whose bounds being received
//! \param[in] lower Lower index of the global rows on my PE
//! \param[in] upper Upper index of the global rows on my PE
// *****************************************************************************
{
  Assert( lower < upper, "Lower bound must be lower than the upper bound: "
          "(" + std::to_string(lower) + "..." +  std::to_string(upper) +
          ") sent by PE " + std::to_string(p) );

  // Store our bounds
  if (p == CkMyPe()) {
    m_lower = lower;
    m_upper = upper;
  }

  // Store inverse of PE-division map stored on all PE
  m_div[ {lower,upper} ] = p;

  // If we have all PEs' bounds, signal the runtime system to continue
  if (m_div.size() == static_cast<std::size_t>(CkNumPes())) {
    // Create my PE's lhs matrix distributed across all PEs
    m_A.create( m_lower*m_ncomp, m_upper*m_ncomp );
    // Create my PE's rhs and unknown vectors distributed across all PEs
    m_b.create( m_lower*m_ncomp, m_upper*m_ncomp );
    m_x.create( m_lower*m_ncomp, m_upper*m_ncomp );
    // Create linear solver
    m_solver.create();
    bounds_complete();
  }
}

void
Solver::next()
// *****************************************************************************
//  Prepare for next step
//! \details Re-enable SDAG waits for rebuilding the right-hand side vector only
// *****************************************************************************
{
  wait4rhsbc();
  wait4hyprerhs();
  wait4asm();
  wait4low();

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

  // Continue with next time step
  for (auto i : m_myworker) m_worker[i].dt();
}

void
Solver::nchare( int n )
// *****************************************************************************
//  Set number of worker chares expected to contribute on my PE
//! \param[in] n Total number of chares (work units) across all PEs
// *****************************************************************************
{
  auto chunksize = n / CkNumPes();
  auto mynchare = chunksize;
  if (CkMyPe() == CkNumPes()-1) mynchare += n % CkNumPes();
  m_nchare = static_cast< std::size_t >( mynchare );

  nchare_complete();
}

void
Solver::charecom( const inciter::CProxy_MatCG& worker,
                  int fromch,
                  const std::vector< std::size_t >& row )
// *****************************************************************************
//  Chares contribute their global row ids for establishing communications
//! \param[in] worker Charm chare array proxy contribution coming from
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] row Global mesh point (row) indices contributed
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
// *****************************************************************************
{
  // Store worker proxy for solution update later
  m_worker = worker;

  // Collect chare ids of workers on my PE
  m_myworker.push_back( fromch );

  // Store rows owned and pack those to be exported, also build import map
  // used to test for completion
  std::map< int, std::set< std::size_t > > exp;
  for (auto gid : row) {
    if (gid >= m_lower && gid < m_upper) {  // if own
      m_rowimport[ fromch ].push_back( gid );
      m_row.insert( gid );
    } else exp[ pe(gid) ].insert( gid );
  }

  // Export non-owned parts to fellow branches that own them
  m_nperow += exp.size();
  for (const auto& p : exp) {
    auto tope = static_cast< int >( p.first );
    thisProxy[ tope ].addrow( fromch, CkMyPe(), p.second );
  }

  if (comcomplete()) contribute( m_cb.get< tag::com >() );
}

void
Solver::addrow( int fromch, int frompe, const std::set< std::size_t >& row )
// *****************************************************************************
//  Receive global row ids from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] frompe PE contribution coming from
//! \param[in] row Global mesh point (row) indices received
// *****************************************************************************
{
  for (auto r : row) {
    m_rowimport[ fromch ].push_back( r );
    m_row.insert( r );
  }
  thisProxy[ frompe ].recrow();
}

void
Solver::recrow()
// *****************************************************************************
//  Acknowledge received row ids
// *****************************************************************************
{
  --m_nperow;
  if (comcomplete()) contribute( m_cb.get< tag::com >() );
}

void
Solver::created()
// *****************************************************************************
// Signal the runtime system that the workers have been created
// *****************************************************************************
{
  if (++m_ncomm == m_nchare) {
    m_ncomm = 0;
    contribute( m_cb.get< tag::com >() );
  }
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
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
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
      exp[ pe(gid[i]) ][ gid[i] ] = solution[i];
    }

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tope = static_cast< int >( p.first );
    thisProxy[ tope ].addsol( fromch, p.second );
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
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
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
      auto& row = exp[ pe(gid[i]) ][ gid[i] ];
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
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
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
      exp[ pe(gid[i]) ][ gid[i] ] = r[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tope = static_cast< int >( p.first );
    thisProxy[ tope ].addrhs( fromch, p.second );
  }

  if (rhscomplete()) rhs_complete();
}

void
Solver::addrhs( int fromch,
                const std::map< std::size_t, std::vector< tk::real > >& r )//,
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
//  Chares contribute to the rhs of the low-order linear system
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] lowrhs Portion of the low-order rhs vector contributed
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
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
      exp[ pe(gid[i]) ][ gid[i] ] = lowrhs[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tope = static_cast< int >( p.first );
    thisProxy[ tope ].addlowrhs( fromch, p.second );
  }

  if (lowrhscomplete()) lowrhs_complete();
}

void
Solver::addlowrhs( int fromch,
                   const std::map< std::size_t,
                                   std::vector< tk::real > >& lowrhs )
// *****************************************************************************
//  Receive+add low-order rhs vector nonzeros from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] lowrhs Portion of the low-order rhs vector contributed,
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
//  Chares contribute to the lhs of the low-order linear system
//! \param[in] fromch Charm chare array the contribution coming from
//! \param[in] gid Global row indices of the vector contributed
//! \param[in] lowlhs Portion of the low-order lhs vector contributed
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
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
      exp[ pe(gid[i]) ][ gid[i] ] = lowlhs[i];

  // Export non-owned vector values to fellow branches that own them
  for (const auto& p : exp) {
    auto tope = static_cast< int >( p.first );
    thisProxy[ tope ].addlowlhs( fromch, p.second );
  }

  if (lowlhscomplete()) lowlhs_complete();
}

void
Solver::addlowlhs( int fromch,
                   const std::map< std::size_t,
                                   std::vector< tk::real > >& lowlhs )
// *****************************************************************************
// Receive and add lhs vector to the low-order system from fellow group branches
//! \param[in] fromch Charm chare array index contribution coming from
//! \param[in] lowlhs Portion of the lhs vector contributed to the low-order
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
//  All communications have been establised among PEs
//! \details At this point all solver objects on all PEs must have received
//!   their global row ids which means that the communications have been
//!   established among all PEs and this the communications (maps) are final on
//!   all PEs.
// *****************************************************************************
{
  // Assert that all global row indices have been received on my PE.  The assert
  // consists of three necessary conditions, which together comprise the
  // sufficient condition that all global row indices have been received owned
  // by this PE.
  Assert( // 1. have heard from every chare on my PE
          m_myworker.size() == m_nchare &&
          // 2. number of rows equals that of the expected on my PE
          m_row.size() == m_upper-m_lower &&
          // 3. all fellow PEs have received my row ids contribution
          m_nperow == 0,
          // if any of the above is unsatisfied, the row ids are incomplete
          "Row ids are incomplete on PE " + std::to_string(CkMyPe()) );

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
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
// *****************************************************************************
{
  // Associate BC vectors to mesh nodes owned
  if (!bc.empty())  // only if chare has anything to offer
    for (const auto& n : bc) {
      Assert( n.second.size() == m_ncomp, "The total number of scalar "
      "components does not equal that of set in the BC vector." );
      m_bc[ n.first ] = n.second;
    }

  // Forward all BC vectors received to fellow branches
  if (++m_nchbc == m_nchare) {
    auto stream = tk::serialize( m_bc );
    m_shadow.ckLocalBranch()->contribute(
      stream.first, stream.second.get(), BCMapMerger,
      CkCallback(CkIndex_Solver::addbc(nullptr),thisProxy) );
  }
}

void
Solver::addbc( CkReductionMsg* msg )
// *****************************************************************************
// Reduction target collecting the final aggregated BC node list map
// *****************************************************************************
{
  PUP::fromMem creator( msg->getData() );
  creator | m_bc;
  delete msg;
  m_nchbc = 0;
  bc_complete();  bc_complete();
}

int
Solver::pe( std::size_t gid )
// *****************************************************************************
//  Return processing element for global mesh row id
//! \param[in] gid Global mesh point (matrix or vector row) id
//! \details First we attempt to the point index in the cache. If that
//!   fails, we resort to a linear search across the division map. Once the
//!   PE is found, we store it in the cache, so next time the search is
//!   quicker. This procedure must find the PE for the id.
//! \return PE that owns global row id
// *****************************************************************************
{
  int p = -1;
  auto it = m_pe.find( gid );
  if (it != end(m_pe))
    p = it->second;
  else
    for (const auto& d : m_div)
      if (gid >= d.first.first && gid < d.first.second)
        p = m_pe[ gid ] = d.second;

  Assert( p >= 0, "PE not found for node id " + std::to_string(gid) );
  Assert( p < CkNumPes(), "Assigning to nonexistent PE" );

  return p;
}

bool
Solver::comcomplete() const
// *****************************************************************************
//  Check if we have done our part in storing and exporting global row ids
//! \details This does not mean the global row ids on our PE is complete
//!   (which is tested by an assert in comfinal), only that we have done
//!   our part of receiving contributions from chare array groups storing
//!   the parts that we own and have sent the parts we do not own to fellow
//!   PEs, i.e., we have nothing else to export. Only when all other fellow
//!   branches have received all contributions are the row ids complete on
//!   all PEs. This latter condition can only be tested after the global
//!   reduction initiated by signal2host_row_complete, which is called when
//!   all fellow branches have returned true from comcomplete.
//! \see comfinal()
//! \return True if we have done our part storing and exporting row ids
// *****************************************************************************
{
  return // have heard from every chare on my PE
         m_myworker.size() == m_nchare &&
         // all fellow PEs have received my row ids contribution
         m_nperow == 0;
}

void
Solver::hyprerow()
// *****************************************************************************
//  Build Hypre data for our portion of the global row ids
//! \note Hypre only likes one-based indexing. Zero-based row indexing fails
//!   to update the vector with HYPRE_IJVectorGetValues().
// *****************************************************************************
{
//std::cout << CkMyPe() << ": " << "hyprerow\n";
  if (m_hypreRows.empty()) {
    for (auto r : m_row) {
      std::vector< int > h( m_ncomp );
      std::iota( begin(h), end(h), r*m_ncomp+1 );
      m_hypreRows.insert( end(m_hypreRows), begin(h), end(h) );
    }
    hyprerow_complete();  hyprerow_complete();  hyprerow_complete();
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
          "Nonzero values of distributed matrix on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );

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
  Assert( rhscomplete(),
          "Values of distributed right-hand-side vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot set BCs" );

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
//  std::cout << CkMyPe() << ": " << "hypresol\n";
  Assert( solcomplete(),
          "Values of distributed solution vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete" );

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
  Assert( lhscomplete(),
          "Nonzero values of distributed matrix on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot convert" );

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
  Assert( rhscomplete(),
          "Values of distributed right-hand-side vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot convert" );

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
  Assert( m_hypreSol.size() == m_hypreRows.size(),
          "Solution vector values incomplete on PE " +
          std::to_string(CkMyPe()) );

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
  Assert( m_hypreMat.size() == m_hypreCols.size(),
          "Matrix values incomplete on PE " + std::to_string(CkMyPe()) );

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
  Assert( m_hypreRhs.size() == m_hypreRows.size(),
          "RHS vector values incomplete on PE " + std::to_string(CkMyPe()) );

  // Set our portion of the vector values
  m_b.set( static_cast< int >( (m_upper - m_lower)*m_ncomp ),
           m_hypreRows.data(),
           m_hypreRhs.data() );

  m_b.assemble();
  asmrhs_complete();
}

void
Solver::updateSol()
// *****************************************************************************
//  Update solution vector in our PE's workers
// *****************************************************************************
{
  // Get solution vector values for our PE
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

    m_worker[ w.first ].updateSol( gid, solution );
  }
}

void
Solver::solve()
// *****************************************************************************
//  Solve hyigh-order linear system
// *****************************************************************************
{
  m_solver.solve( m_A, m_b, m_x );
  updateSol();
}

void
Solver::updateLowSol()
// *****************************************************************************
//  Update low-order solution vector in our PE's workers
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

    m_worker[ w.first ].updateLowSol( gid, solution );
  }
}

void
Solver::lowsolve()
// *****************************************************************************
//  Solve low-order linear system
// *****************************************************************************
{
  // Set boundary conditions on the low-order system
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

  // Solve low-order system
  Assert( rhscomplete(),
          "Values of distributed right-hand-side vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( lowrhscomplete(),
          "Values of distributed mass diffusion rhs vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( lowlhscomplete(),
          "Values of distributed lumped mass lhs vector on PE " +
          std::to_string( CkMyPe() ) + " is incomplete: cannot solve low "
          "order system" );
  Assert( tk::keyEqual( m_rhs, m_lowrhs ), "Row IDs of rhs and mass "
          "diffusion rhs vector unequal on PE " + std::to_string( CkMyPe() )
          + ": cannot solve low order system" );
  Assert( tk::keyEqual( m_rhs, m_lowlhs ), "Row IDs of rhs and lumped mass "
          "lhs vector unequal on PE " + std::to_string( CkMyPe() ) + ": "
          "cannot solve low order system" );

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
