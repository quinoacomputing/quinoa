// *****************************************************************************
/*!
  \file      src/Inciter/DiagCG.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DiagCG for a PDE system with continuous Galerkin without a matrix
  \details   DiagCG advances a system of partial differential equations (PDEs)
    using continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a time
    stepping scheme that is equivalent to the Lax-Wendroff (LW) scheme within
    the unstructured-mesh FE context and treats discontinuities with
    flux-corrected transport (FCT). Only the diagonal entries of the left-hand
    side matrix are non-zero thus it does not need a mtrix-based linear solver.
  \see The documentation in DiagCG.h.
*/
// *****************************************************************************

#include "QuinoaConfig.h"
#include "DiagCG.h"
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
#include "CGPDE.h"
#include "Discretization.h"
#include "DistFCT.h"
#include "DiagReducer.h"
#include "BoundaryConditions.h"

#ifdef HAS_ROOT
  #include "RootMeshWriter.h"
#endif

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< CGPDE > g_cgpde;

} // inciter::

using inciter::DiagCG;

DiagCG::DiagCG( const CProxy_Discretization& disc,
                const tk::CProxy_Solver& solver,
                const FaceData& ) :
  m_itf( 0 ),
  m_nsol( 0 ),
  m_nlhs( 0 ),
  m_nrhs( 0 ),
  m_ndif( 0 ),
  m_disc( disc ),
  m_side( Disc()->BC()->sideNodes( Disc()->Filenodes(), Disc()->Lid() ) ),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_ul( m_u.nunk(), m_u.nprop() ),
  m_du( m_u.nunk(), m_u.nprop() ),
  m_dul( m_u.nunk(), m_u.nprop() ),
  m_ue( m_disc[thisIndex].ckLocal()->Inpoel().size()/4, m_u.nprop() ),
  m_lhs( m_u.nunk(), m_u.nprop() ),
  m_rhs( m_u.nunk(), m_u.nprop() ),
  m_dif( m_u.nunk(), m_u.nprop() ),
  m_bc(),
  m_lhsc(),
  m_rhsc(),
  m_difc(),
  m_vol( 0.0 ),
  m_diag()
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] solver Linear system solver (Solver) proxy
// *****************************************************************************
{
  auto d = Disc();

  // Allocate communication buffers for LHS, ICs, RHS, mass diffusion RHS
  auto np = m_u.nprop();
  auto nb = d->Bid().size();
  m_lhsc.resize( nb );
  for (auto& b : m_lhsc) b.resize( np );
  m_rhsc.resize( nb );
  for (auto& b : m_rhsc) b.resize( np );
  m_difc.resize( nb );
  for (auto& b : m_difc) b.resize( np );

  // Zero communication buffers for setup (LHS, ICs)
  for (auto& b : m_lhsc) std::fill( begin(b), end(b), 0.0 );

  // Signal the runtime system that the workers have been created
  solver.ckLocalBranch()->created();
}

void
DiagCG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types initiated from this chare array
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  NodeDiagnostics::registerReducers();
}

void
DiagCG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  auto d = Disc();

  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Output fields metadata to output file
  d->writeNodeMeta();

  // Compute left-hand side of PDEs
  lhs();

  // Set initial conditions for all PDEs
  for (const auto& eq : g_cgpde) eq.initialize( d->Coord(), m_u, d->T() );

  // Activate SDAG waits for setup
  wait4setup();

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );
}

void
DiagCG::start()
// *****************************************************************************
//  Start time stepping
// *****************************************************************************
{
  auto d = Disc();

  // Start timer measuring time stepping wall clock time
  d->Timer().zero();

  // Combine own and communicated contributions to LHS and ICs
  for (const auto& b : d->Bid()) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    const auto& blhsc = m_lhsc[ b.second ];
    for (ncomp_t c=0; c<m_lhs.nprop(); ++c) m_lhs(lid,c,0) += blhsc[c];
  }

  // Zero communication buffers for first time step (rhs, mass diffusion rhs)
  for (auto& b : m_rhsc) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_difc) std::fill( begin(b), end(b), 0.0 );

  dt();
}

void
DiagCG::lhs()
// *****************************************************************************
// Compute the left-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute lumped mass lhs required for both high and low order solutions
  m_lhs = d->FCT()->lump( *d );

  if (d->Msum().empty())
    comlhs_complete();
  else // send contributions of lhs to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > l;
      for (auto i : n.second) l.push_back( m_lhs[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comlhs( n.second, l );
    }

  ownlhs_complete();
}

void
DiagCG::comlhs( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& L )
// *****************************************************************************
//  Receive contributions to left-hand side diagonal matrix on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive LHS contributions
//! \param[in] L Partial contributions of LHS to chare-boundary nodes
//! \details This function receives contributions to m_lhs, which stores the
//!   diagonal (lumped) mass matrix at mesh nodes. While m_lhs stores
//!   own contributions, m_lhsc collects the neighbor chare contributions during
//!   communication. This way work on m_lhs and m_lhsc is overlapped. The two
//!   are combined in start().
// *****************************************************************************
{
  Assert( L.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  auto d = Disc();

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_lhsc.size(), "Indexing out of bounds" );
    m_lhsc[ bid ] += L[i];
  }

  if (++m_nlhs == d->Msum().size()) {
    m_nlhs = 0;
    comlhs_complete();
  }
}

void
DiagCG::dt()
// *****************************************************************************
// Comppute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  auto d = Disc();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    for (const auto& eq : g_cgpde) {
      auto eqdt = eq.dt( d->Coord(), d->Inpoel(), m_u );
      if (eqdt < mindt) mindt = eqdt;
    }

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Contribute to minimum dt across all chares the advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DiagCG,advance), thisProxy) );
}

void
DiagCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute right-hand side and query Dirichlet BCs for all equations
  for (const auto& eq : g_cgpde)
    eq.rhs( d->T(), d->Dt(), d->Coord(), d->Inpoel(), m_u, m_ue, m_rhs );

  // Query and match user-specified boundary conditions to side sets
  bc();

  if (d->Msum().empty())
    comrhs_complete();
  else // send contributions of rhs to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > r;
      for (auto i : n.second) r.push_back( m_rhs[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comrhs( n.second, r );
    }

  ownrhs_complete();

  // Compute mass diffusion rhs contribution required for the low order solution
  m_dif = d->FCT()->diff( *d, m_u );

  if (d->Msum().empty())
    comdif_complete();
  else // send contributions of diff to chare-boundary nodes to fellow chares
    for (const auto& n : d->Msum()) {
      std::vector< std::vector< tk::real > > D;
      for (auto i : n.second) D.push_back( m_dif[ tk::cref_find(d->Lid(),i) ] );
      thisProxy[ n.first ].comdif( n.second, D );
    }

  owndif_complete();
}

void
DiagCG::comrhs( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& R )
// *****************************************************************************
//  Receive contributions to right-hand side vector on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive RHS contributions
//! \param[in] R Partial contributions of RHS to chare-boundary nodes
//! \details This function receives contributions to m_rhs, which stores the
//!   right hand side vector at mesh nodes. While m_rhs stores own
//!   contributions, m_rhsc collects the neighbor chare contributions during
//!   communication. This way work on m_rhs and m_rhsc is overlapped. The two
//!   are combined in solve().
// *****************************************************************************
{
  Assert( R.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  auto d = Disc();

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_rhsc.size(), "Indexing out of bounds" );
    m_rhsc[ bid ] += R[i];
  }

  if (++m_nrhs == d->Msum().size()) {
    m_nrhs = 0;
    comrhs_complete();
  }
}

void
DiagCG::comdif( const std::vector< std::size_t >& gid,
                const std::vector< std::vector< tk::real > >& D )
// *****************************************************************************
//  Receive contributions to right-hand side mass diffusion on chare-boundaries
//! \param[in] gid Global mesh node IDs at which we receive ontributions
//! \param[in] D Partial contributions to chare-boundary nodes
//! \details This function receives contributions to m_dif, which stores the
//!   mass diffusion right hand side vector at mesh nodes. While m_dif stores
//!   own contributions, m_difc collects the neighbor chare contributions during
//!   communication. This way work on m_dif and m_difc is overlapped. The two
//!   are combined in solve().
// *****************************************************************************
{
  Assert( D.size() == gid.size(), "Size mismatch" );

  using tk::operator+=;

  auto d = Disc();

  for (std::size_t i=0; i<gid.size(); ++i) {
    auto bid = tk::cref_find( d->Bid(), gid[i] );
    Assert( bid < m_difc.size(), "Indexing out of bounds" );
    m_difc[ bid ] += D[i];
  }

  if (++m_ndif == d->Msum().size()) {
    m_ndif = 0;
    comdif_complete();
  }
}

void
DiagCG::bc()
// *****************************************************************************
// Query and match user-specified boundary conditions to side sets
// *****************************************************************************
{
  auto d = Disc();

  // Match user-specified boundary conditions to side sets
  m_bc = d->BC()->match( m_u.nprop(), d->T(), d->Dt(), d->Coord(), d->Gid(),
                         m_side );
}

void
DiagCG::solve()
// *****************************************************************************
//  Solve low and high order diagonal systems
// *****************************************************************************
{
  const auto ncomp = m_rhs.nprop();

  auto d = Disc();

  // Combine own and communicated contributions to rhs and mass diffusion
  for (const auto& b : d->Bid()) {
    auto lid = tk::cref_find( d->Lid(), b.first );
    const auto& brhsc = m_rhsc[ b.second ];
    for (ncomp_t c=0; c<ncomp; ++c) m_rhs(lid,c,0) += brhsc[c];
    const auto& bdifc = m_difc[ b.second ];
    for (ncomp_t c=0; c<ncomp; ++c) m_dif(lid,c,0) += bdifc[c];
  }

  // Zero communication buffers for next time step (rhs, mass diffusion rhs)
  for (auto& b : m_rhsc) std::fill( begin(b), end(b), 0.0 );
  for (auto& b : m_difc) std::fill( begin(b), end(b), 0.0 );

  // Set Dirichlet BCs for lhs and both low and high order rhs vectors. Note
  // that the low order rhs (more prcisely the mass-diffusion term) is set to
  // zero instead of the solution increment at Dirichlet BCs, because for the
  // low order solution the right hand side is the sum of the high order right
  // hand side and mass diffusion so the low order system is L = R + D, where L
  // is the lumped mass matrix, R is the high order RHS, and D is
  // mass diffusion, and R already will have the Dirichlet BC set.
  for (const auto& n : m_bc) {
    auto b = d->Lid().find( n.first );
    if (b != end(d->Lid())) {
      auto id = b->second;
      for (ncomp_t c=0; c<ncomp; ++c)
        if (n.second[c].first) {
          m_lhs( id, c, 0 ) = 1.0;
          m_rhs( id, c, 0 ) = n.second[c].second;
          m_dif( id, c, 0 ) = 0.0;
        }
    }
  }

  // Solve low and high order diagonal systems and update low order solution
  m_dul = (m_rhs + m_dif) / m_lhs;
  m_ul = m_u + m_dul;
  m_du = m_rhs / m_lhs;

  // Continue with FCT
  d->FCT()->aec( *d, m_du, m_u, m_bc );
  d->FCT()->alw( m_u, m_ul, m_dul, thisProxy );  
}

bool
DiagCG::correctBC( const tk::Fields& a )
// *****************************************************************************
//  Verify that the change in the solution at those nodes where Dirichlet
//  boundary conditions are set is exactly the amount the BCs prescribe
//! \param[in] a Limited antidiffusive element contributions
//! \return True if the solution is correct at Dirichlet boundary condition
//!   nodes
// *****************************************************************************
{
  if (m_bc.empty()) return true;

  auto d = Disc();

  // We loop through the map that associates a vector of local node IDs to side
  // set IDs for all side sets read from mesh file. Then for each side set for
  // all mesh nodes on a given side set we attempt to find the global node ID in
  // m_bc, which stores only those nodes (and BC settings) at which the user has
  // configured Dirichlet BCs to be set. Then for all scalar components of all
  // system of systems of PDEs integrated if a BC is to be set for a given
  // component, we compute the low order solution increment + the anti-diffusive
  // element contributions, which is the current solution increment (to be used
  // to update the solution at time n) at that node. This solution increment
  // must equal the BC prescribed at the given node as we solve for solution
  // increments. If not, the BCs are not set correctly, which is an error.
  for (const auto& s : m_side)
    for (auto i : s.second) {
      auto u = m_bc.find( d->Gid()[i] );
      if (u != end(m_bc)) {
        const auto& b = u->second;
        Assert( b.size() == m_u.nprop(), "Size mismatch" );
        for (std::size_t c=0; c<b.size(); ++c)
          if ( b[c].first &&
               std::abs( m_dul(i,c,0) + a(i,c,0) - b[c].second ) >
                 std::numeric_limits< tk::real >::epsilon() ) {
             return false;
          }
      }
    }

  return true;
}

void
DiagCG::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
// *****************************************************************************
{
  auto d = Disc();

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
    for (const auto& eq : g_cgpde) {
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
    d->writeNodeSolution( rmw, m_itf, nodefields() );

  } else
  #endif
  {

    // Create ExodusII writer
    tk::ExodusIIMeshWriter ew( d->OutFilename(), tk::ExoWriter::OPEN );
    // Write time stamp
    ew.writeTimeStamp( m_itf, time );
    // Write node fields to file
    d->writeNodeSolution( ew, m_itf, nodefields() );

  }
}

void
DiagCG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

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
  if ( (std::fabs(d->T() + d->Dt()-term) < eps || (d->It()+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
  {
    writeFields( d->T()+d->Dt() );
  }
}

void
DiagCG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  // Activate SDAG-waits for FCT
  d->FCT()->next();

  // Actiavate SDAG waits for time step
  wait4rhs();

  // Compute rhs for next time step
  rhs();
}

void
DiagCG::next( const tk::Fields& a )
// *****************************************************************************
// Prepare for next step
//! \param[in] a Limited antidiffusive element contributions
// *****************************************************************************
{
  // Apply limited antidiffusive element contributions to low order solution
  if (g_inputdeck.get< tag::discr, tag::fct >())
    m_u = m_ul + a;
  else
    m_u = m_u + m_du;

  auto d = Disc();

  // Output field data to file
  out();
  // Compute diagnostics, e.g., residuals
  auto diag =  m_diag.compute( *d, m_u );
  // Increase number of iterations and physical time
  d->next();
  // Output one-liner status report
  d->status();

  // Evaluate whether to continue with next step
  if (!diag) eval();
}

void
DiagCG::eval()
// *****************************************************************************
// Evaluate whether to continue with next step
// *****************************************************************************
{
  auto d = Disc();

  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();

  // If neither max iterations nor max time reached, continue, otherwise finish
  if (std::fabs(d->T()-term) > eps && d->It() < nstep)
    dt();
  else
    contribute( CkCallback( CkReductionTarget(Transporter,finish), d->Tr() ) );
}

#include "NoWarning/diagcg.def.h"
