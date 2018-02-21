// *****************************************************************************
/*!
  \file      src/Inciter/DG.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.
  \see The documentation in DG.h.
*/
// *****************************************************************************

#include "DG.h"
#include "Discretization.h"
#include "DGPDE.h"
#include "Solver.h"
#include "DiagReducer.h"
#include "Diagnostics.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "ExodusIIMeshWriter.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< DGPDE > g_dgpde;

static CkReduction::reducerType DiagMerger;

} // inciter::

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const tk::CProxy_Solver& solver,
        const FaceData& fd ) :
  m_itf( 0 ),
  m_disc( disc ),
  m_fd( fd ),
  m_nelem( m_disc[thisIndex].ckLocal()->Inpoel().size()/4 ),
  m_u( m_nelem,
       g_inputdeck.get< tag::component >().nprop() ),
  m_vol( 0.0 ),
  m_lhs( m_nelem, 0.0 ),
  m_rhs( m_nelem,
         g_inputdeck.get< tag::component >().nprop() )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
  // Signal the runtime system that the workers have been created
  solver.ckLocalBranch()->created();

  auto d = Disc();

  // Compute face geometry
  m_geoFace = tk::genGeoFaceTri(fd.Ntfac(), fd.Inpofa(), d->Coord());

  // Compute element geometry
  m_geoElem = tk::genGeoElemTet(d->Inpoel(), d->Coord());
}

void
DG::registerReducers()
// *****************************************************************************
//  Configure Charm++ reduction types
//! \details Since this is a [nodeinit] routine, the runtime system executes the
//!   routine exactly once on every logical node early on in the Charm++ init
//!   sequence. Must be static as it is called without an object. See also:
//!   Section "Initializations at Program Startup" at in the Charm++ manual
//!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
// *****************************************************************************
{
  DiagMerger = CkReduction::addReducer( tk::mergeDiag );
}

void
DG::setup( tk::real v )
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
  d->writeElemMeta();

  // Compute left-hand side of discrete PDEs
  lhs();

  // zero initial solution vector
  // ...

  // Set initial conditions for all PDEs
  m_ax = 1.0;
  m_ay = 1.0;
  m_az = 0.0;

  for (std::size_t e=0; e<m_nelem; ++e)
  {
    auto xcc = m_geoElem(e,1,0);
    auto ycc = m_geoElem(e,2,0);
    auto zcc = m_geoElem(e,3,0);
    IGNORE(zcc);

    tk::real u = 1.0 * exp( -((xcc-0.25)*(xcc-0.25) 
                            + (ycc-0.25)*(ycc-0.25))/(2.0 * 0.005) );
    m_u(e,0,0) = u; // scalar
  }

  // Set initial conditions for all PDEs
  std::cout << "DG::init: " << g_dgpde.size() << '\n';
  for (const auto& eq : g_dgpde) eq.initialize( d->Coord(), m_u, d->T() );

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() ) writeFields( d->T() );

  // Start time stepping
  dt();
}

void
DG::dt()
// *****************************************************************************
// Compute time step size
// *****************************************************************************
{
  tk::real mindt = std::numeric_limits< tk::real >::max();

  auto const_dt = g_inputdeck.get< tag::discr, tag::dt >();
  auto def_const_dt = g_inputdeck_defaults.get< tag::discr, tag::dt >();
  auto eps = std::numeric_limits< tk::real >::epsilon();

  // use constant dt if configured
  if (std::abs(const_dt - def_const_dt) > eps) {

    mindt = const_dt;

  } else {      // compute dt based on CFL

    // find the minimum dt across all PDEs integrated
    // ...
    mindt = 0.1;        // stub for now to overwrite numeric_limits::max

    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  // Contribute to minimum dt across all chares the advance to next step
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(DG,advance), thisProxy) );
}

void
DG::writeFields( tk::real time )
// *****************************************************************************
// Output mesh-based fields to file
//! \param[in] time Physical time
// *****************************************************************************
{
  auto d = Disc();

  // Save time stamp at which the last field write happened
  d->LastFieldWriteTime() = time;

  // Increase field output iteration count
  ++m_itf;

  // Collect element field output
  std::vector< std::vector< tk::real > > elemfields;
  elemfields.push_back( m_u.extract(0,0) );

  // Create ExodusII writer
  tk::ExodusIIMeshWriter ew( d->OutFilename(), tk::ExoWriter::OPEN );
  // Write time stamp
  ew.writeTimeStamp( m_itf, time );
  // Write node fields to file
  d->writeElemSolution( ew, m_itf, elemfields );
}

bool
DG::diagnostics()
// *****************************************************************************
// Compute diagnostics, e.g., residuals
//! \return True if diagnostics have been computed
// *****************************************************************************
{
  auto d = Disc();

  const auto ncomp = g_inputdeck.get< tag::component >().nprop();

  std::vector< std::vector< tk::real > >
    diag( NUMDIAG, std::vector< tk::real >( ncomp, 0.0 ) );

  // Compute diagnostics
  // ...

  // Contribute to diagnostics across all PEs
  auto stream = tk::serialize( diag );
  contribute( stream.first, stream.second.get(), DiagMerger,
    CkCallback(CkIndex_Transporter::diagnostics(nullptr), d->Tr()) );

  return true;
}

void
DG::out()
// *****************************************************************************
// Output mesh field data
// *****************************************************************************
{
  auto d = Disc();

  // Optionally output field and particle data
  if ( !((d->It()+1) % g_inputdeck.get< tag::interval, tag::field >()) &&
       !g_inputdeck.get< tag::cmd, tag::benchmark >() )
  {
    writeFields( d->T() + d->Dt() );
  }

  // Output final field data to file (regardless of whether it was requested)
  const auto term = g_inputdeck.get< tag::discr, tag::term >();
  const auto eps = std::numeric_limits< tk::real >::epsilon();
  const auto nstep = g_inputdeck.get< tag::discr, tag::nstep >();
  if ( (std::fabs(d->T() + d->Dt() - term) < eps || (d->It()+1) >= nstep) &&
       (!g_inputdeck.get< tag::cmd, tag::benchmark >()) )
  {
    writeFields( d->T()+d->Dt() );
  }
}

void
DG::lhs()
// *****************************************************************************
// Compute left-hand side of discrete transport equations
// *****************************************************************************
{
  for (std::size_t e=0; e<m_nelem; ++e)
  {
    m_lhs[e] = m_geoElem(e,0,0);
  }
}

void
DG::rhs()
// *****************************************************************************
// Compute right-hand side of discrete transport equations
// *****************************************************************************
{
  auto& esuf = m_fd.Esuf();
  auto& belem = m_fd.Belem();
  auto& triinp = m_fd.Inpofa();

  // initialize rhs as zero
  for (std::size_t e=0; e<m_nelem; ++e)
  {
    m_rhs(e,0,0) = 0.0;
  }

  // compute internal surface flux integrals
  for (auto f=m_fd.Nbfac(); f<m_fd.Ntfac(); ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    auto farea = m_geoFace(f,0,0);

    std::vector< tk::real > fn { m_geoFace(f,1,0),
                                 m_geoFace(f,2,0),
                                 m_geoFace(f,3,0) };

    // need to use tk::Fields::extract() here somehow
    std::vector< tk::real > ul { m_u(el,0,0) };
    std::vector< tk::real > ur { m_u(er,0,0) };

    //--- upwind fluxes
    auto flux = upwindFlux(ul, ur, fn);
    
    m_rhs(el,0,0) -= farea * flux[0];
    m_rhs(er,0,0) += farea * flux[0];
  }

  // compute boundary surface flux integrals
  for (std::size_t f=0; f<m_fd.Nbfac(); ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);

    auto farea = m_geoFace(f,0,0);

    //CkPrintf("(%d/%d) : %d  || ", f, m_fd.Nbfac(), el);
    //CkPrintf("triinpo : %d,  %d,  %d \n", triinp[3*f], triinp[3*f+1], triinp[3*f+2] );

    std::vector< tk::real > fn { m_geoFace(f,1,0),
                                 m_geoFace(f,2,0),
                                 m_geoFace(f,3,0) };

    // need to use tk::Fields::extract() here somehow
    std::vector< tk::real > ul { m_u(el,0,0) };
    std::vector< tk::real > ur { 0.0 };

    //--- upwind fluxes
    auto flux = upwindFlux(ul, ur, fn);
    
    m_rhs(el,0,0) -= farea * flux[0];
  }
}

std::vector< tk::real >
DG::upwindFlux( std::vector< tk::real > ul,
                std::vector< tk::real > ur,
                std::vector< tk::real > fn )
// *****************************************************************************
// Riemann solver using upwind method
//! \param[in] ul Left unknown/state vector
//! \param[in] ur Right unknown/state vector
//! \param[in] fn Face unit normal vector
//! \return Riemann solution using upwind method
// *****************************************************************************
{
    std::vector< tk::real > flux(ul.size(),0);

    // wave speed
    tk::real swave = m_ax*fn[0] + m_ay*fn[1] + m_az*fn[2];

    // upwinding
    if (swave > 0.0)
    {
      flux[0] = swave * ul[0];
    }
    else
    {
      flux[0] = swave * ur[0];
    }

    return flux;
}

void
DG::tstep()
// *****************************************************************************
// Explicit time-stepping using forward Euler to discretize time-derivative
// *****************************************************************************
{
}

void
DG::advance( tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
// *****************************************************************************
{
  auto d = Disc();

  // Set new time step size
  d->setdt( newdt );

  // Compute rhs for next time step
  rhs();

  // Advance solution/time-stepping
  tstep();

  // Prepare for next time step
  next();
}

void
DG::next()
// *****************************************************************************
// Prepare for next step
// *****************************************************************************
{
  auto d = Disc();

  // Output field data to file
  out();
  // Compute diagnostics, e.g., residuals
  auto diag = diagnostics();
  // Increase number of iterations and physical time
  d->next();
  // Output one-liner status report
  d->status();

  // Evaluate whether to continue with next step
  if (!diag) eval();
}

void
DG::eval()
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

#include "NoWarning/dg.def.h"
