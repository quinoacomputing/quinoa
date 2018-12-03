// *****************************************************************************
/*!
  \file      src/PDE/FunctionPrototypes.h
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Function prototypes used in PDE classes
  \details   This file defines prototypes used in PDE classes. Examples include
     functions used for evaluating known (e.g., analytical) solutions or setting
     initial conditions (ICs), and flux functions used for evaluating Riemann
     fluxes for discontinuous Galerkin methods.
*/
// *****************************************************************************
#ifndef FunctionPrototypes_h
#define FunctionPrototypes_h

#include <vector>
#include <functional>

#include "Types.h"
#include "Keywords.h"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Function prototype for Problem::solution() functions
//! \details Functions of this type are used to evaluate known (e.g.,
//!    analytical) solutions or setting initial conditions
//! \see e.g., inciter::CompFlowProblemVorticalFlow::solution
//! \note Used for both contininuous and discontinuous Galerkin discretizations
using SolutionFn = std::function<
  std::vector< real >( ncomp_t, ncomp_t, real, real, real, real ) >;

//! Function prototype for Riemann flux functions
//! \details Functions of this type are used to compute numerical fluxes across a
//!    surface
//! \see e.g., inciter::Upwind, inciter::LaxFriedrichs, inciter::HLLC
using FluxFn = std::function<
  std::vector< real >( const std::array< real, 3 >&,
                       const std::array< std::vector< real >, 2 >&,
                       const std::vector< std::array< real, 3 > >& ) >;

//! Function prototype for evaluating a prescribed velocity field
//! \details Functions of this type are used to prescribe known velocity fields
//! \note Used for scalar transport
//! \see e.g., TransportProblemShearDiff::prescribedVelocity
using VelFn = std::function<
  std::vector< std::array< tk::real, 3 > >
  ( real, real, real, ncomp_t, ncomp_t ) >;

} // tk::

#endif // FunctionPrototypes_h
