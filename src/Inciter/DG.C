// *****************************************************************************
/*!
  \file      src/Inciter/DG.C
  \copyright 2012-2015, J. Bakosi, 2016-2017, Los Alamos National Security, LLC.
  \brief     DG advances a system of PDEs with the discontinuous Galerkin scheme
  \details   DG advances a system of partial differential equations (PDEs) using
    discontinuous Galerkin (DG) finite element (FE) spatial discretization (on
    tetrahedron elements) combined with Runge-Kutta (RK) time stepping.
  \see The documentation in DG.h.
*/
// *****************************************************************************

#include "DG.h"
#include "Discretization.h"
#include "PDE.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_inputdeck_defaults;
extern std::vector< PDE > g_pdes;

} // inciter::

using inciter::DG;

DG::DG( const CProxy_Discretization& disc,
        const tk::CProxy_Solver& ) :
  m_disc( disc ),
  m_vol( 0.0 )
// *****************************************************************************
//  Constructor
// *****************************************************************************
{
}

void
DG::setup( tk::real v )
// *****************************************************************************
// Setup rows, query boundary conditions, output mesh, etc.
//! \param[in] v Total mesh volume
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();
  // Store total mesh volume
  m_vol = v;
  // Output chare mesh to file
  d->writeMesh();
  // Output fields metadata to output file
  d->writeMeta();
}

void
DG::init()
// *****************************************************************************
// Set ICs, compute initial time step size, output initial field data
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();

  // zero initial solution vector
  // ...

  // Set initial conditions for all PDEs
  // ...

  // Compute initial time step size
  dt();

  // Output initial conditions to file (regardless of whether it was requested)
  if ( !g_inputdeck.get< tag::cmd, tag::benchmark >() )
    writeFields( d->T() );
}

void
DG::writeFields( tk::real )
// *****************************************************************************
// Output mesh-based fields to file
// *****************************************************************************
{
}

void
DG::dt()
// *****************************************************************************
// Comppute time step size
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
    /// ...
 
    // Scale smallest dt with CFL coefficient
    mindt *= g_inputdeck.get< tag::discr, tag::cfl >();

  }

  auto d = m_disc[ thisIndex ].ckLocal();

  // Contribute to mindt across all CG chares
  contribute( sizeof(tk::real), &mindt, CkReduction::min_double,
              CkCallback(CkReductionTarget(Transporter,dt), d->Tr()) );
}

void
DG::advance( uint64_t it, tk::real t, tk::real newdt )
// *****************************************************************************
// Advance equations to next time step
//! \param[in] newdt Size of this new time step
//! \param[in] it Iteration count
//! \param[in] t Physical time
// *****************************************************************************
{
  auto d = m_disc[ thisIndex ].ckLocal();

  // Update local copy of time step info (the master copies are in Transporter)
  d->It() = it;
  d->T() = t;
  d->Dt() = newdt;

  // Compute rhs for next time step
  // ...
}

#include "NoWarning/dg.def.h"
