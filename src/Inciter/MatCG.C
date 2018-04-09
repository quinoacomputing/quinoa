// *****************************************************************************
/*!
  \file      src/Inciter/MatCG.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     MatCG for a PDE system with continuous Galerkin with a matrix
  \details   MatCG advances a system of partial differential equations (PDEs)
    using continuous Galerkin (CG) finite element (FE) spatial discretization
    (using linear shapefunctions on tetrahedron elements) combined with a time
    stepping scheme that is equivalent to the Lax-Wendroff (LW) scheme within
    the unstructured-mesh FE context and treats discontinuities with
    flux-corrected transport (FCT). The left-hand side matrix is stored in a
    compressed sparse row (CSR) storage and thus uses a matrix-based linear
    solver.
  \see The documentation in MatCG.h.
*/
// *****************************************************************************

#include "QuinoaConfig.h"
#include "MatCG.h"
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

using inciter::MatCG;

MatCG::MatCG( const CProxy_Discretization& disc,
              const tk::CProxy_Solver& solver,
              const FaceData& ) :
  m_itf( 0 ),
  m_nhsol( 0 ),
  m_nlsol( 0 ),
  m_disc( disc ),
  m_solver( solver ),
  m_side( Disc()->BC()->sideNodes( Disc()->Filenodes(), Disc()->Lid() ) ),
  m_u( m_disc[thisIndex].ckLocal()->Gid().size(),
       g_inputdeck.get< tag::component >().nprop() ),
  m_ul( m_u.nunk(), m_u.nprop() ),
  m_du( m_u.nunk(), m_u.nprop() ),
  m_dul( m_u.nunk(), m_u.nprop() ),
  m_ue( m_disc[thisIndex].ckLocal()->Inpoel().size()/4, m_u.nprop() ),
  m_lhsd( m_disc[thisIndex].ckLocal()->Psup().second.size()-1, m_u.nprop() ),
  m_lhso( m_disc[thisIndex].ckLocal()->Psup().first.size(), m_u.nprop() ),
  m_vol( 0.0 ),
  m_diag( *Disc() )
// *****************************************************************************
//  Constructor
//! \param[in] disc Discretization proxy
//! \param[in] solver Linear system solver (Solver) proxy
// *****************************************************************************
{
  auto d = Disc();

  // Send off global row IDs to linear system solver
  m_solver.ckLocalBranch()->charecom( thisProxy, thisIndex, d->Gid() );
}

void
MatCG::registerReducers()
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
MatCG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  auto d = Disc();

  m_solver.ckLocalBranch()->comfinal();

  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Output fields metadata to output file
  d->writeNodeMeta();

  // Compute left-hand side of PDEs
  lhs();

  // zero initial solution vector
  m_du.fill( 0.0 );

  // Set initial conditions for all PDEs
  for (const auto& eq : g_cgpde) eq.initialize( d->Coord(), m_u, d->T() );

  // Send off initial guess for assembly
  m_solver.ckLocalBranch()->charesol( thisIndex, d->Gid(), m_du );

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );

  // Start time stepping
  contribute( CkCallback( CkReductionTarget(Transporter,start), d->Tr()) );

  // Start timer measuring time stepping wall clock time
  d->Timer().zero();
}

void
MatCG::dt()
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
              CkCallback(CkReductionTarget(MatCG,advance), thisProxy) );
}

void
MatCG::lhs()
// *****************************************************************************
// Compute left-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute left-hand side matrix for all equations
  for (const auto& eq : g_cgpde)
    eq.lhs( d->Coord(), d->Inpoel(), d->Psup(), m_lhsd, m_lhso );

  // Send off left hand side for assembly
  m_solver.ckLocalBranch()->
    charelhs( thisIndex, d->Gid(), d->Psup(), m_lhsd, m_lhso );

  // Compute lumped mass lhs required for the low order solution
  auto lump = d->FCT()->lump( *d );
  // Send off lumped mass lhs for assembly
  m_solver.ckLocalBranch()->charelowlhs( thisIndex, d->Gid(), lump );
}

void
MatCG::rhs()
// *****************************************************************************
// Compute right-hand side of transport equations
// *****************************************************************************
{
  auto d = Disc();

  // Compute right-hand side and query Dirichlet BCs for all equations
  tk::Fields r( d->Gid().size(), g_inputdeck.get< tag::component >().nprop() );
  for (const auto& eq : g_cgpde)
    eq.rhs( d->T(), d->Dt(), d->Coord(), d->Inpoel(), m_u, m_ue, r );

  // Query and match user-specified boundary conditions to side sets
  bc();

  // Send off right-hand sides for assembly
  m_solver.ckLocalBranch()->charerhs( thisIndex, d->Gid(), r );

  // Compute mass diffusion rhs contribution required for the low order solution
  auto diff = d->FCT()->diff( *d, m_u );
  // Send off mass diffusion rhs contribution for assembly
  m_solver.ckLocalBranch()->charelowrhs( thisIndex, d->Gid(), diff );
}

void
MatCG::bc()
// *****************************************************************************
// Query and match user-specified boundary conditions to side sets
// *****************************************************************************
{
  auto d = Disc();

  // Match user-specified boundary conditions to side sets
  auto dirbc = d->BC()->match( m_u.nprop(), d->T(), d->Dt(), d->Coord(),
                               d->Gid(), m_side );

  // Send off BCs to Solver for aggregation
  m_solver.ckLocalBranch()->charebc( dirbc );
}

void
MatCG::updateLowSol( const std::vector< std::size_t >& gid,
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

  auto d = Disc();

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( d->Lid(), gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_dul( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nlsol += gid.size();

  // If all contributions we own have been received, continue with FCT
  if (m_nlsol == d->Gid().size()) {
    m_nlsol = 0;
    m_ul = m_u + m_dul;
    d->FCT()->alw( m_u, m_ul, m_dul, thisProxy );
  }
}

void
MatCG::updateSol( const std::vector< std::size_t >& gid,
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

  auto d = Disc();

  // Receive update of solution vector
  for (std::size_t i=0; i<gid.size(); ++i) {
    auto id = tk::cref_find( d->Lid(), gid[i] );
    for (ncomp_t c=0; c<ncomp; ++c) m_du( id, c, 0 ) = du[ i*ncomp+c ];
  }

  // Count number of solution nodes updated
  m_nhsol += gid.size();

  // If all contributions we own have been received, continue with FCT
  if (m_nhsol == d->Gid().size()) {
    m_nhsol = 0;
    d->FCT()->aec( *d, m_du, m_u, m_solver.ckLocalBranch()->dirbc() );
  }
}

bool
MatCG::correctBC( const tk::Fields& a )
// *****************************************************************************
//  Verify that the change in the solution at those nodes where Dirichlet
//  boundary conditions are set is exactly the amount the BCs prescribe
//! \param[in] a Limited antidiffusive element contributions
//! \return True if the solution is correct at Dirichlet boundary condition
//!   nodes
// *****************************************************************************
{
  auto& dirbc = m_solver.ckLocalBranch()->dirbc();

  if (dirbc.empty()) return true;

  auto d = Disc();

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
               std::abs( m_dul(i,c,0) + a(i,c,0) - b[c].second ) >
                 std::numeric_limits< tk::real >::epsilon() ) {
             return false;
          }
      }
  }

  return true;
}

void
MatCG::writeFields( tk::real time )
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
MatCG::out()
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
MatCG::advance( tk::real newdt )
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

  // Compute rhs for next time step
  rhs();
}

void
MatCG::next( const tk::Fields& a )
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
  auto diag = m_diag.compute( *d, m_u );
  // Increase number of iterations and physical time
  d->next();
  // Output one-liner status report
  d->status();

//     // TEST FEATURE: Manually migrate this chare by using migrateMe to see if
//     // all relevant state variables are being PUPed correctly.
//     //CkPrintf("I'm MatCG chare %d on PE %d\n",thisIndex,CkMyPe());
//     if (thisIndex == 2 && CkMyPe() == 2) {
//       /*int j;
//       for (int i; i < 50*std::pow(thisIndex,4); i++) {
//         j = i*thisIndex;
//       }*/
//       CkPrintf("I'm MatCG chare %d on PE %d\n",thisIndex,CkMyPe());
//       migrateMe(1);
//    }
//    if (thisIndex == 2 && CkMyPe() == 1) {
//      CkPrintf("I'm MatCG chare %d on PE %d\n",thisIndex,CkMyPe());
//      migrateMe(2);
//    }

  // Evaluate whether to continue with next step
  if (!diag) eval();
}

void
MatCG::eval()
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
    contribute( CkCallback( CkReductionTarget(Transporter,next), d->Tr() ) );
  else
    contribute( CkCallback( CkReductionTarget(Transporter,finish), d->Tr() ) );
}

#include "NoWarning/matcg.def.h"
