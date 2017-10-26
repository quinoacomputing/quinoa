// *****************************************************************************
/*!
  \file      src/Inciter/CG.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     CG advances a system of PDEs with the continuous Galerkin scheme
  \details   CG advances a system of partial differential equations (PDEs) using
    continuous Galerkin (CG) finite element (FE) spatial discretization (using
    linear shapefunctions on tetrahedron elements) combined with a time stepping
    scheme that is equivalent to the Lax-Wendroff (LW) scheme within the
    unstructured-mesh FE context and treats discontinuities with flux-corrected
    transport (FCT).
  \see The documentation in CG.h.
*/
// *****************************************************************************

#include <string>
#include <cmath>
#include <array>
#include <set>
#include <algorithm>

#include "QuinoaConfig.h"
#include "CG.h"
#include "Solver.h"
#include "Vector.h"
#include "Reader.h"
#include "ContainerUtil.h"
#include "UnsMesh.h"
#include "Reorder.h"
#include "ExodusIIMeshReader.h"
#include "ExodusIIMeshWriter.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "DerivedData.h"
#include "PDE.h"
#include "Discretization.h"

#ifdef HAS_ROOT
  #include "RootMeshWriter.h"
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< PDE > g_pdes;

} // inciter::

using inciter::CG;

CG::CG( const CProxy_Discretization& disc,
        const tk::CProxy_Solver& solver ) :
  m_itf( 0 ),
  m_nhsol( 0 ),
  m_nlsol( 0 ),
  m_naec( 0 ),
  m_nalw( 0 ),
  m_nlim( 0 ),
  m_disc( disc ),
  m_solver( solver ),
  m_side(),
  m_fluxcorrector( m_disc[thisIndex].ckLocal()->Inpoel().size() ),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_ul( m_u.nunk(), m_u.nprop() ),
  m_du( m_u.nunk(), m_u.nprop() ),
  m_dul( m_u.nunk(), m_u.nprop() ),
  m_ue( m_disc[thisIndex].ckLocal()->Inpoel().size()/4, m_u.nprop() ),
  m_p( m_u.nunk(), m_u.nprop()*2 ),
  m_q( m_u.nunk(), m_u.nprop()*2 ),
  m_a( m_u.nunk(), m_u.nprop() ),
  m_lhsd( m_disc[thisIndex].ckLocal()->Psup().second.size()-1, m_u.nprop() ),
  m_lhso( m_disc[thisIndex].ckLocal()->Psup().first.size(), m_u.nprop() ),
  m_pc(),
  m_qc(),
  m_ac(),
  m_vol( 0.0 )
// *****************************************************************************
//  Constructor
//! \param[in] transporter Host (Transporter) proxy
//! \param[in] disc Discretization proxy
//! \param[in] solver Linear system solver (Solver) proxy
//! \param[in] filenodes Map associating old node IDs (as in file) to new node
//!   IDs (as in producing contiguous-row-id linear system contributions)
//! \details "Contiguous-row-id" here means that the numbering of the mesh nodes
//!   (which corresponds to rows in the linear system) are (approximately)
//!   contiguous (as much as this can be done with an unstructured mesh) as the
//!   problem is distirbuted across PEs, held by Solver objects. This ordering
//!   is in start contrast with "as-in-file" ordering, which is the ordering of
//!   the mesh nodes as it is stored in the file from which the mesh is read in.
//!   The as-in-file ordering is highly non-contiguous across the distributed
//!   problem.
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Allocate receive buffers for FCT
  m_pc.resize( d->Bid().size() );
  for (auto& b : m_pc) b.resize( m_u.nprop()*2 );
  m_qc.resize( d->Bid().size() );
  for (auto& b : m_qc) b.resize( m_u.nprop()*2 );
  m_ac.resize( d->Bid().size() );
  for (auto& b : m_ac) b.resize( m_u.nprop() );

  // Invert file-node map, a map associating old node IDs (as in file) to new
  // node IDs (as in producing contiguous-row-id linear system contributions),
  // so we can search more efficiently for old node IDs.
  std::unordered_map< std::size_t, std::size_t > linnodes;
  for (const auto& i : d->Filenodes()) {
    auto n = tk::cref_find(d->Lid(),i.first);
    Assert( n < d->Gid().size(),
            "Local IDs must be lower than the local number of grid points" );
    linnodes[ i.second ] = n;
  }

  // Access all side sets and their old node IDs (as in file) from Solver
  auto& oldside = m_solver.ckLocalBranch()->side();

  // Create map that assigns the local mesh node IDs mapped to side set ids,
  // storing only those nodes for a given side set that are part of our chunk of
  // the mesh.
  for (const auto& s : oldside) {
    auto& n = m_side[ s.first ];
    for (auto o : s.second) {
      auto it = linnodes.find( o );
      if (it != end(linnodes))
        n.push_back( it->second );
    }
  }

  // Signal the runtime system that the CG worker objects have been created
  contribute(CkCallback(CkReductionTarget(Transporter,created), d->Tr()));
}

void
CG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );
  
  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Send off global row IDs to linear system solver
  m_solver.ckLocalBranch()->charecom( thisIndex, d->Gid() );
  // Output fields metadata to output file
  d->writeMeta();
}

void
CG::bc()
// *****************************************************************************
//  Extract node IDs from side set node lists and match to user-specified
//  boundary conditions
//! \details Boundary conditions (BC), mathematically speaking, are applied on
//!   finite surfaces. These finite surfaces are given by element sets (i.e., a
//!   list of elements). This function queries Dirichlet boundary condition
//!   values from all PDEs in the system of PDEs integrated at the node lists
//!   associated to side set IDs (previously read from file). As a response to
//!   this query, each PDE system returns a BC data structure which is then
//!   sent to the linear system solver which needs to know about this to apply
//!   BCs before a linear solve. Note that the BC mesh nodes that this function
//!   results in, stored in dirbc and sent to the linear system solver, only
//!   contains those nodes that this chare contributes to, i.e., it does not
//!   contain those BC nodes at which other chares enforce Dirichlet BCs. The
//!   linear system solver then collects these and communicates to other PEs so
//!   that BC data held in Solver::m_bc are the same on all PEs. That is
// *****************************************************************************
{
  using NodeBC = std::vector< std::pair< bool, tk::real > >;

  // Vector of pairs of bool and boundary condition value associated to mesh
  // node IDs at which the user has set Dirichlet boundary conditions for all
  // PDEs integrated
  std::unordered_map< std::size_t, NodeBC > dirbc;

  // Query Dirichlet boundary conditions for all PDEs integrated and assign to
  // nodes. This is where the individual system of PDEs are queried for boundary
  // conditions. The outer loop goes through all sides sets that exists in the
  // input file and passes the map's value_type (a pair of the side set id and a
  // vector of local node IDs) to PDE::dirbc(). PDE::dirbc() returns a new map
  // that associates a vector of pairs associated to local node IDs. (The pair
  // is a pair of bool and real value, the former is the fact that the BC is to
  // be set while the latter is the value if it is to be set). The length of
  // this NodeBC vector, returning from each system of PDEs equals to the number
  // of scalar components the given PDE integrates. Here then we contatenate
  // this map for all PDEs integrated. If there are multiple BCs set at a mesh
  // node (dirbc::key), either because (1) in the same PDE system the user
  // prescribed BCs on side sets that share nodes or (2) because more than a
  // single PDE system assigns BCs to a given node (on different variables), the
  // NodeBC vector must be correctly stored. "Correctly" here means that the
  // size of the NodeBC vectors must all be the same and qual to the sum of all
  // scalar components integrated by all PDE systems integrated. Example:
  // single-phase compressible flow (density, momentum, energy = 5) +
  // transported scalars of 10 variables -> NodeBC vector length = 15. Note that
  // in case (1) above a new node encountered must "overwrite" the already
  // existing space for the NodeBC vector. "Overwrite" here means that it should
  // keep the existing BCs and add the new ones yielding the union the two
  // prescription for BCs but in the same space that already exist in the NodeBC
  // vector. In case (2), however, the NodeBC pairs must go to the location in
  // the vector assigned to the given PDE system, i.e., using the above example
  // BCs for the 10 (or less) scalars should go in the positions starting at 5,
  // leaving the first 5 false, indicating no BCs for the flow variables.
  //
  // TODO: Note that the logic described above is only partially implemented at
  // this point. What works is the correct insertion of multiple BCs for nodes
  // shared among multiple side sets, e.g., corners, originating from the same
  // PDE system. What is not yet implemented is the case when there are no BCs
  // set for flow variables but there are BCs for transport, the else branch
  // below will incorrectly NOT skip the space for the flow variables. In other
  // words, this only works for a single PDE system and a sytem of systems. This
  // machinery is only tested with a single system of PDEs at this point.

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  for (const auto& s : m_side) {
    std::size_t c = 0;
    for (std::size_t eq=0; eq<g_pdes.size(); ++eq) {
      auto eqbc = g_pdes[eq].dirbc( d->T(), d->Dt(), s, d->Coord() );
      for (const auto& n : eqbc) {
        auto id = n.first;                      // BC node ID
        const auto& bcs = n.second;             // BCs
        auto& nodebc = dirbc[ d->Gid()[id] ];   // BCs to be set for node
        if (nodebc.size() > c) {        // node already has BCs from this PDE
          Assert( nodebc.size() == c+bcs.size(), "Size mismatch" );
          for (std::size_t i=0; i<bcs.size(); i++)
            if (bcs[i].first)
              nodebc[c+i] = bcs[i];
        } else {        // node does not yet have BCs from this PDE
          // This branch needs to be completed for system of systems of PDEs.
          // See note above.
          nodebc.insert( end(nodebc), begin(bcs), end(bcs) );
        }
      }
      if (!eqbc.empty()) c += eqbc.cbegin()->second.size();
    }
  }

  // Verify the size of each NodeBC vectors. They must have the same lengths and
  // equal to the total number of scalar components for all systems of PDEs
  // integrated. This is intentional, because this way the linear system solver
  // does not have to (and does not) know about individual equation systems.
  for (const auto& n : dirbc) {
    IGNORE(n);
    Assert( n.second.size() == m_u.nprop(), "Size of NodeBC vector incorrect" );
  }

//   // Send progress report to host
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chbcmatched();

  // Send off list of owned node IDs mapped to side sets to Solver
  m_solver.ckLocalBranch()->charebc( dirbc );
}

void
CG::init()
// *****************************************************************************
// Set ICs, compute initial time step size, output initial field data, compute
// left-hand-side matrix
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // zero initial solution vector
  m_du.fill( 0.0 );

  // Send off initial guess for assembly
  m_solver.ckLocalBranch()->charesol( thisIndex, d->Gid(), m_du );

  // Set initial conditions for all PDEs
  for (const auto& eq : g_pdes)
    eq.initialize( d->Coord(), m_u, d->T(), d->Gid() );

  // Compute initial time step size
  dt();

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );

//   // send progress report to host
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chic();

  // Compute left-hand side of PDEs
  lhs();
}

void
CG::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Only write if the last time is different than the current one
  if (std::abs(d->LastFieldWriteTime() - time) <
      std::numeric_limits< tk::real >::epsilon() )
    return;

  // Save time stamp at which the last field write happened
  d->LastFieldWriteTime() = time;

  // Increase field output iteration count
  ++m_itf;

  // Lambda to collect node fields output from all PDEs
  auto nodefields = [&]() {
    auto u = m_u;   // make a copy as eq::output() may overwrite its arg
    std::vector< std::vector< tk::real > > output;
    for (const auto& eq : g_pdes) {
      auto o = eq.fieldOutput( time, m_vol, d->Coord(), d->V(), u );
      output.insert( end(output), begin(o), end(o) );
    }
    return output;
  };

  #ifdef HAS_ROOT
  auto filetype = g_inputdeck.get< tag::selected, tag::filetype >();

  if (filetype == tk::ctr::FieldFileType::ROOT) {

    // Create Root writer
    tk::RootMeshWriter rmw( d->OutFilename(), 1 );
    // Write time stamp
    rmw.writeTimeStamp( m_itf, time );
    // Write node fields to file
    d->writeSolution( rmw, m_itf, nodefields() );

  } else
  #endif
  {

    // Create ExodusII writer
    tk::ExodusIIMeshWriter ew( d->OutFilename(), tk::ExoWriter::OPEN );
    // Write time stamp
    ew.writeTimeStamp( m_itf, time );
    // Write node fields to file
    d->writeSolution( ew, m_itf, nodefields() );

  }
}

void
CG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Optionally output field and particle data
  if ( !((d->It()+1) % g_inputdeck.get< tag::interval, tag::field >()) &&
       !g_inputdeck.get< tag::cmd, tag::benchmark >() )
  {
    writeFields( d->T()+d->Dt() );
  }

  // Output final field data to file (regardless of whether it was requested)
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  if ( (std::fabs(d->T()+d->Dt()-term) < eps || (d->It()+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
    writeFields( d->T()+d->Dt() );
}


void
CG::dt()
// *****************************************************************************
// Comppute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    for (const auto& eq : g_pdes) {
      auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_u );
      if (eqdt < mindt) mindt = eqdt;
    }

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Contribute to mindt across all CG chares
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(Transporter,dt), d->Tr()) );
}

void
CG::lhs()
// *****************************************************************************
// Compute left-hand side of transport equations
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_pdes)
    eq.lhs( d->Coord(), d->Inpoel(), d->Psup(), m_lhsd, m_lhso );

  // Send off left hand side for assembly
  m_solver.ckLocalBranch()->
    charelhs( thisIndex, d->Gid(), d->Psup(), m_lhsd, m_lhso );

  // Compute lumped mass lhs required for the low order solution
  auto lump = m_fluxcorrector.lump( d->Coord(), d->Inpoel() );
  // Send off lumped mass lhs for assembly
  m_solver.ckLocalBranch()->charelowlhs( thisIndex, d->Gid(), lump );

//   // send progress report to host
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chlhs();
}

void
CG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  // Initialize FCT data structures for new time step
  m_p.fill( 0.0 );
  m_a.fill( 0.0 );
  for (std::size_t p=0; p<m_u.nunk(); ++p)
    for (ncomp_t c=0; c<g_inputdeck.get<tag::component>().nprop(); ++c) {
      m_q(p,c*2+0,0) = -std::numeric_limits< tk::real >::max();
      m_q(p,c*2+1,0) = std::numeric_limits< tk::real >::max();
    }

  for (auto& b : m_pc) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_ac) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_qc)
    for (ncomp_t c=0; c<m_u.nprop(); ++c) {
      b[c*2+0] = -std::numeric_limits< tk::real >::max();
      b[c*2+1] = std::numeric_limits< tk::real >::max();
    }

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Compute right-hand side and query Dirichlet BCs for all equations
  tk::Fields r( d->Gid().size(), g_inputdeck.get< tag::component >().nprop() );
  for (const auto& eq : g_pdes)
    eq.rhs( d->T(), d->Dt(), d->Coord(), d->Inpoel(), m_u, m_ue, r );
  // Query Dirichlet BCs and send to linear system solver
  bc();
  // Send off right-hand sides for assembly
  m_solver.ckLocalBranch()->charerhs( thisProxy, thisIndex, d->Gid(), r );

  // Compute mass diffusion rhs contribution required for the low order solution
  auto diff = m_fluxcorrector.diff( d->Coord(), d->Inpoel(), m_u );
  // Send off mass diffusion rhs contribution for assembly
  m_solver.ckLocalBranch()->charelowrhs( thisIndex, d->Gid(), diff );

//   // send progress report to host
//   auto d = m_disc[ thisIndex ].ckLocal();
//   Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chrhs();
}

void
CG::aec()
// *****************************************************************************
//  Compute and sum antidiffusive element contributions (AEC) to mesh nodes
//! \details This function computes and starts communicating m_p, which stores
//!    the sum of all positive (negative) antidiffusive element contributions to
//!    nodes (Lohner: P^{+,-}_i), see also FluxCorrector::aec().
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Compute and sum antidiffusive element contributions to mesh nodes. Note
  // that the sums are complete on nodes that are not shared with other chares
  // and only partial sums on chare-boundary nodes.
  auto& dbc = m_solver.ckLocalBranch()->dirbc();
  m_fluxcorrector.aec( d->Coord(), d->Inpoel(), d->Vol(), dbc, d->Gid(), m_du,
                       m_u, m_p );

  if (d->Msum().empty())
    comaec_complete();
  else // send contributions to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > p;
      for (auto i : n.second) p.push_back( m_p[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comaec( n.second, p );
    }

  ownaec_complete();
  #ifndef NDEBUG
  ownaec_complete();
  #endif
}

void
CG::comaec( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& P )
// *****************************************************************************
//  Receive sums of antidiffusive element contributions on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive AEC contributions
//! \param[in] P Partial sums of positive (negative) antidiffusive element
//!   contributions to chare-boundary nodes
//! \details This function receives contributions to m_p, which stores the
//!   sum of all positive (negative) antidiffusive element contributions to
//!   nodes (Lohner: P^{+,-}_i), see also FluxCorrector::aec(). While m_p stores
//!   own contributions, m_pc collects the neighbor chare contributions during
//!   communication. This way work on m_p and m_pc is overlapped. The two are
//!   combined in lim().
// *****************************************************************************
{
  Assert( P.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_pc.size(), "Indexing out of bounds" );
    m_pc[ bid ] += P[i];
  }

  if (++m_naec == d->Msum().size()) {
    m_naec = 0;
    comaec_complete();
  }
}

void
CG::alw()
// *****************************************************************************
//  Compute the maximum and minimum unknowns of elements surrounding nodes
//! \details This function computes and starts communicating m_q, which stores
//!    the maximum and mimimum unknowns of all elements surrounding each node
//!    (Lohner: u^{max,min}_i), see also FluxCorrector::alw().
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Compute the maximum and minimum unknowns of all elements surrounding nodes
  // Note that the maximum and minimum unknowns are complete on nodes that are
  // not shared with other chares and only partially complete on chare-boundary
  // nodes.
  m_fluxcorrector.alw( d->Inpoel(), m_u, m_ul, m_q );

  if (d->Msum().empty())
    comalw_complete();
  else // send contributions at chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > q;
      for (auto i : n.second) q.push_back( m_q[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comalw( n.second, q );
    }

  ownalw_complete();
  #ifndef NDEBUG
  ownalw_complete();
  #endif
}

void
CG::comalw( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& Q )
// *****************************************************************************
// Receive contributions to the maxima and minima of unknowns of all elements
// surrounding mesh nodes on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] Q Partial contributions to maximum and minimum unknowns of all
//!   elements surrounding nodes to chare-boundary nodes
//! \details This function receives contributions to m_q, which stores the
//!   maximum and mimimum unknowns of all elements surrounding each node
//!   (Lohner: u^{max,min}_i), see also FluxCorrector::alw(). While m_q stores
//!   own contributions, m_qc collects the neighbor chare contributions during
//!   communication. This way work on m_q and m_qc is overlapped. The two are
//!   combined in lim().
// *****************************************************************************
{
  Assert( Q.size() == gid.size(), "Size mismatch" );

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_qc.size(), "Indexing out of bounds" );
    auto& o = m_qc[ bid ];
    const auto& q = Q[i];
    for (ncomp_t c=0; c<m_u.nprop(); ++c) {
      if (q[c*2+0] > o[c*2+0]) o[c*2+0] = q[c*2+0];
      if (q[c*2+1] < o[c*2+1]) o[c*2+1] = q[c*2+1];
    }
  }

  if (++m_nalw == d->Msum().size()) {
    m_nalw = 0;
    comalw_complete();
  }
}

void
CG::lim()
// *****************************************************************************
//  Compute the limited antidiffusive element contributions
//! \details This function computes and starts communicating m_a, which stores
//!   the limited antidiffusive element contributions assembled to nodes
//!   (Lohner: AEC^c), see also FluxCorrector::limit().
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Combine own and communicated contributions to P and Q
  for (const auto& b : d->Bid()) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    const auto& bpc = m_pc[ b.second ];
    const auto& bqc = m_qc[ b.second ];
    for (ncomp_t c=0; c<m_p.nprop()/2; ++c) {
      m_p(lid,c*2+0,0) += bpc[c*2+0];
      m_p(lid,c*2+1,0) += bpc[c*2+1];
      if (bqc[c*2+0] > m_q(lid,c*2+0,0)) m_q(lid,c*2+0,0) = bqc[c*2+0];
      if (bqc[c*2+1] < m_q(lid,c*2+1,0)) m_q(lid,c*2+1,0) = bqc[c*2+1];
    }
  }

  m_fluxcorrector.lim( d->Inpoel(), m_p, m_ul, m_q, m_a );

  if (d->Msum().empty())
    comlim_complete();
  else // send contributions to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > a;
      for (auto i : n.second) a.push_back( m_a[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comlim( n.second, a );
    }

  ownlim_complete();
}

void
CG::comlim( const std::vector< std::size_t >& gid,
                 const std::vector< std::vector< tk::real > >& A )
// *****************************************************************************
//  Receive contributions of limited antidiffusive element contributions on
//  chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive contributions
//! \param[in] A Partial contributions to antidiffusive element contributions to
//!   chare-boundary nodes
//! \details This function receives contributions to m_a, which stores the
//!   limited antidiffusive element contributions assembled to nodes (Lohner:
//!   AEC^c), see also FluxCorrector::limit(). While m_a stores own
//!   contributions, m_ac collects the neighbor chare contributions during
//!   communication. This way work on m_a and m_ac is overlapped. The two are
//!   combined in apply().
// *****************************************************************************
{
  Assert( A.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_ac.size(), "Indexing out of bounds" );
    m_ac[ bid ] += A[i];
  }
 
  if (++m_nlim == d->Msum().size()) {
    m_nlim = 0;
    comlim_complete();
  }
}

void
CG::advance( uint64_t it, tk::real t, tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Update local copy of time step info (the master copies are in Transporter)
  d->It() = it;
  d->T() = t;
  d->Dt() = newdt;

  // Activate SDAG-waits
  #ifndef NDEBUG
  wait4ver();
  #endif
  wait4fct();
  wait4app();

  // Compute rhs for next time step
  rhs();
}

void
CG::updateLowSol( const std::vector< std::size_t >& gid,
                  const std::vector< tk::real >& du )
// *****************************************************************************
// Update low order solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] du Portion of the unknown/solution vector update
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  Assert( gid.size() * ncomp == du.size(),
          "Size of row ID vector times the number of scalar components and the "
          "size of the low order solution vector must equal" );

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( d->Lid(), gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_dul( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nlsol += gid.size();

  // If all contributions we own have been received, continue by updating a
  // different solution vector
  if (m_nlsol == d->Gid().size()) {
    m_nlsol = 0;
    m_ul = m_u + m_dul;
    alw();
  }
}

void
CG::updateSol( const std::vector< std::size_t >& gid,
               const std::vector< tk::real >& du )
// *****************************************************************************
// Update high order solution vector
//! \param[in] gid Global row indices of the vector updated
//! \param[in] du Portion of the unknown/solution vector update
// *****************************************************************************
{
  auto ncomp = g_inputdeck.get< tag::component >().nprop();
  Assert( gid.size() * ncomp == du.size(),
          "Size of row ID vector times the number of scalar components and the "
          "size of the high order solution vector must equal" );

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( d->Lid(), gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_du( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nhsol += gid.size();

  // If all contributions we own have been received, continue by updating a
  // different solution vector
  if (m_nhsol == d->Gid().size()) {
    m_nhsol = 0;
    aec();
  }
}

void
CG::verify()
// *****************************************************************************
// Verify antidiffusive element contributions up to linear solver convergence
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  if (m_fluxcorrector.verify( d->Nchare(), d->Inpoel(), m_du, m_dul ))
    contribute( CkCallback( CkReductionTarget(Transporter,verified), d->Tr()) );
}

void
CG::diagnostics()
// *****************************************************************************
// Compute diagnostics, e.g., residuals
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Optionally: collect analytical solutions and send both the latest
  // analytical and numerical solutions to Solver for computing and outputing
  // diagnostics
  if ( !(d->It() % g_inputdeck.get< tag::interval, tag::diag >()) ) {

    // Collect analytical solutions (if available) from all PDEs. Note that
    // calling the polymorphic PDE::initialize() is assumed to evaluate the
    // analytical solution for a PDE. For those PDE problems that have
    // analytical solutions, this is the same as used for setting the initial
    // conditions, since if the analytical solution is available for a problem
    // it is (so far anyway) always initialized to its analytical solution and
    // that is done by calling PDE::initialize(). If the analytical solution is
    // a function of time, that is already incorporated in setting initial
    // conditions. For those PDEs where the analytical solution is not
    // available, initialize() returns the initial conditions (obviously), and
    // thus the "error", defined between this "analytical" and numerical
    // solution will be a measure of the "distance" between the initial
    // condition and the current numerical solution. This is not necessarily
    // useful, but simplies the logic because all PDEs can be treated as baing
    // able to compute an error based on some "analytical" solution, which
    // really the initial condition.
    auto a = m_u;
    for (const auto& eq : g_pdes)
      eq.initialize( d->Coord(), a, d->T()+d->Dt(), d->Gid() );

    // Prepare for computing diagnostics. Diagnostics are defined as the L2 norm
    // of a quantity, computed in mesh nodes, A, as || A ||_2 = sqrt[ sum_i (
    // A_i )^2 V_i ], where the sum is taken over all mesh nodes and V_i is the
    // nodal volume. We send two sets of quantities to the host for aggregation
    // across the whole mesh: (1) the numerical solutions of all components of
    // all PDEs, and their error, defined as A_i = (a_i - n_i), where a_i and
    // n_i are the analytical and numerical solutions at node i, respectively.
    // We send these to Solver and the final aggregated solution will end up in
    // Transporter::diagnostics(). Note that we weigh/multiply all data here by
    // sqrt(V_i), so that the nodal volumes do not have to be communicated
    // separately. In Solver::diagnostics(), where we collect all contributions
    // from chares on a PE, all data is squared. Solver::diagnostics() is where
    // the sums are computed, then the sums are summed accross the whole problem
    // in Transporter::diagnostics(), where the final square-root of the L2
    // norm, defined above, is taken.

    // Send both numerical and analytical solutions to solver
    m_solver.ckLocalBranch()->charediag( thisIndex, d->Gid(), m_u, a, d->V() );

  } else
    contribute(
      CkCallback(CkReductionTarget(Transporter,diagcomplete), d->Tr()) );
}

bool
CG::correctBC()
// *****************************************************************************
//  Verify that the change in the solution at those nodes where Dirichlet
//  boundary conditions are set is exactly the amount the BCs prescribe
//! \return True if the solution is correct at Dirichlet boundary condition
//!   nodes
// *****************************************************************************
{
  auto& dirbc = m_solver.ckLocalBranch()->dirbc();

  if (dirbc.empty()) return true;

  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // We loop through the map that associates a vector of local node IDs to side
  // set IDs for all side sets read from mesh file. Then for each side set for
  // all mesh nodes on a given side set we attempt to find the global node ID
  // in dirbc, which stores only those nodes (and BC settings) at which the
  // user has configured Dirichlet BCs to be set. Then for all scalar
  // components of all system of systems of PDEs integrated if a BC is to be
  // set for a given component, we compute the low order solution increment +
  // the anti-diffusive element contributions, which is the current solution
  // increment (to be used to update the solution at time n) at that node. This
  // solution increment must equal the BC prescribed at the given node as we
  // solve for solution increments. If not, the BCs are not set correctly,
  // which is an error.
  for (const auto& s : m_side)
    for (auto i : s.second) {
      auto u = dirbc.find( d->Gid()[i] );
      if (u != end(dirbc)) {
        const auto& b = u->second;
        Assert( b.size() == m_u.nprop(), "Size mismatch" );
        for (std::size_t c=0; c<b.size(); ++c)
          if ( b[c].first &&
               std::abs( m_dul(i,c,0) + m_a(i,c,0) - b[c].second ) >
                 std::numeric_limits< tk::real >::epsilon() ) {
             return false;
          }
      }
  }

  return true;
}

void
CG::apply()
// *****************************************************************************
// Apply limited antidiffusive element contributions
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  Assert( d!=nullptr, "Discretization proxy's ckLocal() null" );

  // Combine own and communicated contributions to A
  for (const auto& b : d->Bid()) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    const auto& bac = m_ac[ b.second ];
    for (ncomp_t c=0; c<m_a.nprop(); ++c) m_a(lid,c,0) += bac[c];
  }

  // Verify that solution values do not change at Dirichlet BC nodes
  Assert( correctBC(), "Dirichlet boundary condition incorrect" );

  // Apply limited antidiffusive element contributions to low order solution
  if (g_inputdeck.get< tag::discr, tag::fct >())
    m_u = m_ul + m_a;
  else
    m_u = m_u + m_du;

//   // send progress report to host
//   if ( g_inputdeck.get< tag::cmd, tag::feedback >() ) d->Tr().chlim();

  // Output field data to file
  out();
  // Compute diagnostics, e.g., residuals
  diagnostics();

//     // TEST FEATURE: Manually migrate this chare by using migrateMe to see if
//     // all relevant state variables are being PUPed correctly.
//     //CkPrintf("I'm CG chare %d on PE %d\n",thisIndex,CkMyPe());
//     if (thisIndex == 2 && CkMyPe() == 2) {
//       /*int j;
//       for (int i; i < 50*std::pow(thisIndex,4); i++) {
//         j = i*thisIndex;
//       }*/
//       CkPrintf("I'm CG chare %d on PE %d\n",thisIndex,CkMyPe());
//       migrateMe(1);
//    }
//    if (thisIndex == 2 && CkMyPe() == 1) {
//      CkPrintf("I'm CG chare %d on PE %d\n",thisIndex,CkMyPe());
//      migrateMe(2);
//    }
}

#include "NoWarning/cg.def.h"
