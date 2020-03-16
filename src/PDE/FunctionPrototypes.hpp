// *****************************************************************************
/*!
  \file      src/PDE/FunctionPrototypes.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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

#include "Types.hpp"
#include "Keywords.hpp"
#include "Fields.hpp"

namespace tk {

using ncomp_t = kw::ncomp::info::expect::type;

//! Function prototype for Problem::solution() functions
//! \details Functions of this type are used to evaluate known (e.g.,
//!    analytical) solutions or setting initial conditions
//! \see e.g., inciter::CompFlowProblemVorticalFlow::solution
//! \note Used for both continuous and discontinuous Galerkin discretizations
using SolutionFn = std::function<
  std::vector< real >( ncomp_t, ncomp_t, real, real, real, real, int& ) >;

//! Function prototype for Riemann flux functions
//! \details Functions of this type are used to compute numerical fluxes across a
//!    surface using a Riemann solver
//! \see e.g., inciter::Upwind, inciter::LaxFriedrichs, inciter::HLLC
using RiemannFluxFn = std::function<
  std::vector< real >( const std::array< real, 3 >&,
                       const std::array< std::vector< real >, 2 >&,
                       const std::vector< std::array< real, 3 > >& ) >;

//! Function prototype for flux vector functions
//! \details Functions of this type are used to compute physical flux functions
//!   in the PDEs being solved. These are different than the RiemannFluxFn
//!   because they compute the actual flux functions, not the solution to a
//!   Riemann problem.
//! \see e.g., inciter::dg::Transport::flux, inciter::dg::CompFlow::flux
using FluxFn = std::function<
  std::vector< std::array< real, 3 > >
  ( ncomp_t, ncomp_t, const std::vector< real >&,
    const std::vector< std::array< real, 3 > >& ) >;

//! Function prototype for evaluating a prescribed velocity field
//! \details Functions of this type are used to prescribe known velocity fields
//! \note Used for scalar transport
//! \see e.g., TransportProblemShearDiff::prescribedVelocity
using VelFn = std::function<
  std::vector< std::array< tk::real, 3 > >
  ( ncomp_t, ncomp_t, real, real, real ) >;

//! Function prototype for physical boundary states
//! \details Functions of this type are used to provide the left and right
//!    states of boundary faces along physical boundaries
using StateFn = std::function<
  std::array< std::vector< real >, 2 >
  ( ncomp_t, ncomp_t, const std::vector< real >&, real, real, real, real,
    const std::array< tk::real, 3 >& ) >;

//! Function prototype for evaluating a source term for a system of PDEs
//! \details Functions of this type are used to evaluate an arbitrary source
//!   term specialized to a particular problem, e.g., derived using the method
//!   of manufactured solutions, and specialized to a particular system of PDEs
//! \see e.g., CompFlowProblemRayleighTaylor::src
using SrcFn = std::function<
  std::vector< tk::real >( ncomp_t, ncomp_t, real, real, real, real ) >;

//! \brief Function prototype for computing the element gradient contribution to a
//!    nodal gradient in ALECG
using ElemGradFn = std::function<
  std::tuple< std::array< std::size_t, 4 >,
              std::array< std::array< tk::real, 3 >, 4 >,
              std::vector< std::array< tk::real, 4 > >,
              tk::real >( ncomp_t,
                          ncomp_t,
                          std::size_t,
                          const std::array< std::vector< tk::real >, 3 >&,
                          const std::vector< std::size_t >&,
                          const std::tuple< std::vector< tk::real >,
                                            std::vector< tk::real > >&,
                          const tk::Fields& ) >;

} // tk::

#endif // FunctionPrototypes_h
