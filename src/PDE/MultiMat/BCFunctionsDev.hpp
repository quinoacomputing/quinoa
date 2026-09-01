// *****************************************************************************
/*!
  \file      src/PDE/MultiMat/BCFunctionsDev.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2023 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Device-callable boundary state functions for the _constP path
  \details   Mirrors of the symmetry and extrapolate state functions in
             BCFunctions.hpp, for use by the const-order boundary surface
             integral. The originals are untouched, so the p-adaptive path is
             unaffected.

             Differences from the host versions, all forced by device
             compilation:
               - the right state is written into a caller-owned fixed-size
                 buffer instead of a returned std::vector (the host version
                 does `auto ur = ul`, a heap allocation per face per
                 quadrature point)
               - nmat, ncomp, nsld and solidx are parameters rather than
                 input-deck lookups
               - std::array   -> Kokkos::Array
               - tk::reflectTensor -> tk::reflectTensorDev (BLAS-free)
               - stressCmp[i][j]   -> inciter::stressCmpDev(i,j)
             The arithmetic and its ordering are otherwise unchanged.
*/
// *****************************************************************************
#ifndef BCFunctionsDev_h
#define BCFunctionsDev_h

#include "Riemann/RiemannDevice.hpp"
#include "MultiMat/MultiMatIndexing.hpp"

namespace inciter {

//! Which boundary state function a _constP boundary integral should apply
//! \details A runtime tag rather than a template parameter: a single
//!   bndSurfInt_constP call handles exactly one boundary condition, so the
//!   branch is uniform across every face and quadrature point in the launch
//!   and costs no divergence.
enum class BCKind : std::uint8_t { Symmetry, Extrapolate, NoSlipWall };

//! Extrapolate boundary state: the ghost state equals the internal state
template< class StateT >
KOKKOS_INLINE_FUNCTION void
extrapolateDev( std::size_t nstate, const StateT& ul, StateT& ur )
{
  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];
}

//! Shared wall-type boundary state body
//! \details symmetry and noslipwall differ only in how the ghost velocity is
//!   formed; everything else, including the solid reflection, is identical in
//!   the host originals. The ghost velocity is passed in so that one body
//!   serves both and cannot drift.
template< class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
wallStateDev( std::size_t ncomp,
              std::size_t nmat,
              std::size_t nstate,
              const SolidxT& solidx,
              const FnT& fn,
              const StateT& ul,
              StateT& ur,
              tk::real rho,
              tk::real v1r,
              tk::real v2r,
              tk::real v3r )
{
  // Boundary condition
  for (std::size_t k=0; k<nmat; ++k)
  {
    ur[volfracIdx(nmat, k)] = ul[volfracIdx(nmat, k)];
    ur[densityIdx(nmat, k)] = ul[densityIdx(nmat, k)];
    ur[energyIdx(nmat, k)] = ul[energyIdx(nmat, k)];
    if (solidx[k] > 0) {
      // Internal inverse deformation tensor
      tk::Mat3Dev g;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          g[i][j] = ul[deformIdx(nmat,solidx[k],i,j)];
      // Internal Cauchy stress tensor
      tk::Mat3Dev s;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          s[i][j] = ul[ncomp+stressIdx(nmat,solidx[k],stressCmpDev(i,j))];
      // Make reflection matrix
      tk::Mat3Dev reflectionMat;
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j)
          reflectionMat[i][j] = (i==j ? 1.0 : 0.0) - 2*fn[i]*fn[j];
      // Reflect g and s
      g = tk::reflectTensorDev(g, reflectionMat);
      s = tk::reflectTensorDev(s, reflectionMat);
      // Copy g and s into ur
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j) {
          ur[deformIdx(nmat,solidx[k],i,j)] = g[i][j];
          ur[ncomp+stressIdx(nmat,solidx[k],stressCmpDev(i,j))] = s[i][j];
        }
    }
  }

  ur[momentumIdx(nmat, 0)] = rho * v1r;
  ur[momentumIdx(nmat, 1)] = rho * v2r;
  ur[momentumIdx(nmat, 2)] = rho * v3r;

  // velocity
  ur[ncomp+velocityIdx(nmat, 0)] = v1r;
  ur[ncomp+velocityIdx(nmat, 1)] = v2r;
  ur[ncomp+velocityIdx(nmat, 2)] = v3r;
  // material pressures
  for (std::size_t k=0; k<nmat; ++k)
    ur[ncomp+pressureIdx(nmat, k)] = ul[ncomp+pressureIdx(nmat, k)];
}

//! Symmetry boundary state: velocity reflected about the face normal
template< class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
symmetryDev( std::size_t ncomp,
             std::size_t nmat,
             std::size_t nstate,
             const SolidxT& solidx,
             const FnT& fn,
             const StateT& ul,
             StateT& ur )
{
  tk::real rho(0.0);
  for (std::size_t k=0; k<nmat; ++k)
    rho += ul[densityIdx(nmat, k)];

  // host does `auto ur = ul`; copy explicitly into the caller's buffer
  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];

  // Internal cell velocity components
  auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
  auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
  auto v3l = ul[ncomp+velocityIdx(nmat, 2)];
  // Normal component of velocity
  auto vnl = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];
  // Ghost state velocity components
  auto v1r = v1l - 2.0*vnl*fn[0];
  auto v2r = v2l - 2.0*vnl*fn[1];
  auto v3r = v3l - 2.0*vnl*fn[2];

  wallStateDev( ncomp, nmat, nstate, solidx, fn, ul, ur, rho, v1r, v2r, v3r );
}

//! No-slip wall boundary state: velocity negated
//! \details Identical to symmetry apart from the ghost velocity, matching the
//!   host noslipwall in BCFunctions.hpp.
template< class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
noslipwallDev( std::size_t ncomp,
               std::size_t nmat,
               std::size_t nstate,
               const SolidxT& solidx,
               const FnT& fn,
               const StateT& ul,
               StateT& ur )
{
  tk::real rho(0.0);
  for (std::size_t k=0; k<nmat; ++k)
    rho += ul[densityIdx(nmat, k)];

  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];

  // Internal cell velocity components
  auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
  auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
  auto v3l = ul[ncomp+velocityIdx(nmat, 2)];
  // Ghost state velocity components
  auto v1r = -v1l;
  auto v2r = -v2l;
  auto v3r = -v3l;

  wallStateDev( ncomp, nmat, nstate, solidx, fn, ul, ur, rho, v1r, v2r, v3r );
}

//! Dispatch to the boundary state function selected by kind
template< class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
bcStateDev( BCKind kind,
            std::size_t ncomp,
            std::size_t nmat,
            std::size_t nstate,
            const SolidxT& solidx,
            const FnT& fn,
            const StateT& ul,
            StateT& ur )
{
  if (kind == BCKind::Symmetry)
    symmetryDev( ncomp, nmat, nstate, solidx, fn, ul, ur );
  else if (kind == BCKind::NoSlipWall)
    noslipwallDev( ncomp, nmat, nstate, solidx, fn, ul, ur );
  else
    extrapolateDev( nstate, ul, ur );
}

} // inciter::

#endif // BCFunctionsDev_h
