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
#include "EoS/EOS.hpp"

namespace inciter {

//! Which boundary state function a _constP boundary integral should apply
//! \details A runtime tag rather than a template parameter: a single
//!   bndSurfInt_constP call handles exactly one boundary condition, so the
//!   branch is uniform across every face and quadrature point in the launch
//!   and costs no divergence.
enum class BCKind : std::uint8_t { Symmetry, Extrapolate, NoSlipWall,
                                  Inlet, Farfield, BackPressure };

//! \brief Scalar boundary condition parameters, read from the input deck
//! \details The input deck is not reachable from device code, so the scalars
//!   that inlet, farfield and back_pressure look up on the host are gathered
//!   here and passed into the kernel. A POD struct rather than a growing
//!   parameter list, so adding a boundary condition does not change
//!   bndSurfIntMultiMat_constP's signature again.
//! \note The host inlet reads g_inputdeck bc[0].inlet[0] and the host farfield
//!   reads bc[0], i.e. neither supports per-sideset parameters. That
//!   limitation is carried over here unchanged rather than generalised.
struct BCParamsDev {
  //! multimat min_volumefrac
  tk::real alphamin{0.0};
  //! farfield pressure, temperature, velocity and material index
  tk::real fp{0.0};
  tk::real ft{0.0};
  tk::real fu[3]{0.0, 0.0, 0.0};
  std::size_t fmat{0};
  //! inlet pressure, temperature, velocity and material index
  tk::real p_in{0.0};
  tk::real t_in{0.0};
  tk::real u_in[3]{0.0, 0.0, 0.0};
  std::size_t mat_in{0};
  //! back pressure
  tk::real fbp{0.0};
};

//! Whether a boundary state function needs EOS::density / EOS::totalenergy
//! \details Symmetry, Extrapolate and NoSlipWall derive the ghost state
//!   entirely from the internal state, so they need neither. Inlet, Farfield
//!   and BackPressure all call compute< EOS::density > and
//!   compute< EOS::totalenergy >; BackPressure additionally calls
//!   compute< EOS::temperature >.
//!
//!   This deliberately returns true for any BCKind not listed below, so a
//!   newly added kind is treated as needing density until someone states
//!   otherwise. That is the fail-safe direction: the cost of a wrong "true" is
//!   that the boundary condition runs on the host, whereas the cost of a wrong
//!   "false" is a device kernel silently producing NaN for a material whose
//!   density() is host-only (currently JWL). Do not add a default case.
KOKKOS_INLINE_FUNCTION bool bcKindNeedsDensity( BCKind kind )
{
  switch (kind) {
    case BCKind::Symmetry:
    case BCKind::Extrapolate:
    case BCKind::NoSlipWall:
      return false;
    case BCKind::Inlet:
    case BCKind::Farfield:
    case BCKind::BackPressure:
      return true;
  }
  return true;
}

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

//! Inlet boundary state: velocity, pressure and material from the deck
//! \details Mirror of inciter::inlet in BCFunctions.hpp. g_inputdeck lookups
//!   become BCParamsDev fields; std::max -> fmax.
//! \note The host passes no deformation gradient to compute< EOS::totalenergy >,
//!   letting it default to an empty std::array. The device overload requires an
//!   explicit argument, so a zero matrix is passed to match. For StiffenedGas
//!   that is identical (defgrad is ignored). For a solid material the host's
//!   empty matrix flows into elasticEnergy -> getIsochorRightCauchyGreen and
//!   yields NaN; that is a pre-existing host behaviour, reproduced here rather
//!   than corrected, so host and device agree.
template< class MatBlkT, class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
inletDev( std::size_t ncomp,
          std::size_t nmat,
          std::size_t nstate,
          const SolidxT& /*solidx*/,
          const FnT& fn,
          const StateT& ul,
          StateT& ur,
          const MatBlkT& mat_blk,
          const BCParamsDev& bp )
{
  // host does `auto ur = ul`; copy explicitly into the caller's buffer
  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];

  // External cell velocity, such that velocity = v_in at face
  auto v1r = bp.u_in[0];
  auto v2r = bp.u_in[1];
  auto v3r = bp.u_in[2];

  // Normal inlet velocity
  auto vn = bp.u_in[0]*fn[0] + bp.u_in[1]*fn[1] + bp.u_in[2]*fn[2];

  // Acoustic speed
  tk::real a(0.0);
  for (std::size_t k=0; k<nmat; ++k)
    if (ul[volfracIdx(nmat, k)] > 1.0e-04)
      a = fmax( a, mat_blk[k].template compute< EOS::soundspeed >(
        ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
        ul[volfracIdx(nmat, k)], k ) );

  // Mach number
  auto Ma = vn / a;

  // see the note above: stands in for the host's defaulted empty defgrad
  const tk::real zerog[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

  tk::real pk(0.0);
  tk::real rho(0.0);
  for (std::size_t k=0; k<nmat; ++k) {
    if (k == bp.mat_in)
      ur[volfracIdx(nmat,k)] = 1.0 -
        (static_cast< tk::real >(nmat-1))*bp.alphamin;
    else
      ur[volfracIdx(nmat,k)] = bp.alphamin;

    // Material pressure, which, for supersonic inflow, is the exterior
    // pressure and the interior pressure for subsonic
    if(Ma <= -1)
      pk = bp.p_in;
    else
      pk = ul[ncomp+pressureIdx(nmat,k)]/ul[volfracIdx(nmat,k)];
    auto rhok = mat_blk[k].template compute< EOS::density >(pk, bp.t_in);

    ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat,k)] * pk;
    ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;
    ur[energyIdx(nmat,k)] =
      mat_blk[k].template compute< EOS::totalenergy >(
      ur[volfracIdx(nmat,k)]*rhok,
      v1r, v2r, v3r, ur[volfracIdx(nmat,k)]*pk, ur[volfracIdx(nmat,k)],
      zerog);

    // bulk density
    rho += ur[densityIdx(nmat,k)];
  }

  ur[momentumIdx(nmat, 0)] = rho * v1r;
  ur[momentumIdx(nmat, 1)] = rho * v2r;
  ur[momentumIdx(nmat, 2)] = rho * v3r;

  // velocity
  ur[ncomp+velocityIdx(nmat, 0)] = v1r;
  ur[ncomp+velocityIdx(nmat, 1)] = v2r;
  ur[ncomp+velocityIdx(nmat, 2)] = v3r;
}

//! Farfield boundary state, from the characteristic theory of hyperbolic systems
//! \details Mirror of inciter::farfield in BCFunctions.hpp. All four Mach
//!   regimes are reproduced, including the supersonic-outflow case which the
//!   host handles by leaving ur equal to ul. g_inputdeck lookups become
//!   BCParamsDev fields; std::max -> fmax; tk::getDeformGrad ->
//!   getDeformGradDev; stressCmp[i][j] -> stressCmpDev(i,j).
//! \warning The subsonic-outflow and both inflow branches call
//!   compute< EOS::totalenergy > with a deformation gradient, so for solid
//!   materials this routes through elasticEnergyDev and is NOT bit-identical
//!   to the host. See the warning in EoS/TensorEOSDev.hpp.
template< class MatBlkT, class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
farfieldDev( std::size_t ncomp,
             std::size_t nmat,
             std::size_t nstate,
             const SolidxT& solidx,
             const FnT& fn,
             const StateT& ul,
             StateT& ur,
             const MatBlkT& mat_blk,
             const BCParamsDev& bp )
{
  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];

  // Internal cell velocity components
  auto v1l = ul[ncomp+velocityIdx(nmat, 0)];
  auto v2l = ul[ncomp+velocityIdx(nmat, 1)];
  auto v3l = ul[ncomp+velocityIdx(nmat, 2)];

  // Normal velocity
  auto vn = v1l*fn[0] + v2l*fn[1] + v3l*fn[2];

  // Acoustic speed
  tk::real a(0.0);
  for (std::size_t k=0; k<nmat; ++k)
    if (ul[volfracIdx(nmat, k)] > 1.0e-04)
      a = fmax( a, mat_blk[k].template compute< EOS::soundspeed >(
        ul[densityIdx(nmat, k)], ul[ncomp+pressureIdx(nmat, k)],
        ul[volfracIdx(nmat, k)], k ) );

  // Mach number
  auto Ma = vn / a;

  // Inverse deformation gradient
  tk::real gk[3][3];

  if (Ma <= -1) {  // Supersonic inflow
    // For supersonic inflow, all the characteristics are from outside.
    // Therefore, we calculate the ghost cell state using the primitive
    // variables from outside.
    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      if (k == bp.fmat)
        ur[volfracIdx(nmat,k)] = 1.0 -
          (static_cast< tk::real >(nmat-1))*bp.alphamin;
      else
        ur[volfracIdx(nmat,k)] = bp.alphamin;
      auto rhok = mat_blk[k].template compute< EOS::density >(bp.fp, bp.ft);
      ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;

      // solids' state
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) {
            // inverse deformation gradient
            if (i == j) gk[i][j] = 1.0;
            else gk[i][j] = 0.0;
            ur[deformIdx(nmat,solidx[k],i,j)] = gk[i][j];

            // elastic component of Cauchy stress
            ur[ncomp+stressIdx(nmat,solidx[k],stressCmpDev(i,j))] = 0.0;
          }
      }
      else {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) gk[i][j] = 0.0;
      }

      ur[energyIdx(nmat,k)] =
        mat_blk[k].template compute< EOS::totalenergy >(
        ur[volfracIdx(nmat,k)]*rhok,
        bp.fu[0], bp.fu[1], bp.fu[2], ur[volfracIdx(nmat,k)]*bp.fp,
        ur[volfracIdx(nmat,k)], gk);

      // material pressures
      ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat, k)] * bp.fp;

      rho += ur[densityIdx(nmat,k)];
    }
    for (std::size_t i=0; i<3; ++i) {
      ur[momentumIdx(nmat,i)] = rho * bp.fu[i];
      ur[ncomp+velocityIdx(nmat, i)] = bp.fu[i];
    }

  } else if (Ma > -1 && Ma < 0) {  // Subsonic inflow
    // For subsonic inflow, there is 1 outgoing characteristic and 4
    // incoming characteristics. Therefore, we calculate the ghost cell state
    // by taking pressure from the internal cell and other quantities from
    // the outside.
    tk::real rho(0.0);
    for (std::size_t k=0; k<nmat; ++k) {
      if (k == bp.fmat)
        ur[volfracIdx(nmat,k)] = 1.0 -
          (static_cast< tk::real >(nmat-1))*bp.alphamin;
      else
        ur[volfracIdx(nmat,k)] = bp.alphamin;
      auto p = ul[ncomp+pressureIdx(nmat,k)] / ul[volfracIdx(nmat,k)];
      auto rhok = mat_blk[k].template compute< EOS::density >(p, bp.ft);
      ur[densityIdx(nmat,k)] = ur[volfracIdx(nmat,k)] * rhok;

      // solids' state
      if (solidx[k] > 0) {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) {
            // inverse deformation gradient
            if (i == j) gk[i][j] = 1.0;
            else gk[i][j] = 0.0;
            ur[deformIdx(nmat,solidx[k],i,j)] = gk[i][j];

            // elastic component of Cauchy stress
            ur[ncomp+stressIdx(nmat,solidx[k],stressCmpDev(i,j))] = 0.0;
          }
      }
      else {
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j) gk[i][j] = 0.0;
      }

      ur[energyIdx(nmat,k)] =
        mat_blk[k].template compute< EOS::totalenergy >(
        ur[volfracIdx(nmat,k)]*rhok,
        bp.fu[0], bp.fu[1], bp.fu[2], ur[volfracIdx(nmat,k)]*p,
        ur[volfracIdx(nmat,k)], gk);

      // material pressures
      ur[ncomp+pressureIdx(nmat, k)] = ur[volfracIdx(nmat, k)] * p;

      rho += ur[densityIdx(nmat,k)];
    }
    for (std::size_t i=0; i<3; ++i) {
      ur[momentumIdx(nmat,i)] = rho * bp.fu[i];
      ur[ncomp+velocityIdx(nmat, i)] = bp.fu[i];
    }

  } else if (Ma >= 0 && Ma < 1) {  // Subsonic outflow
    // For subsonic outflow, there is 1 incoming characteristic and 4
    // outgoing characteristics. Therefore, we calculate the ghost cell state
    // by taking pressure from the outside and other quantities from the
    // internal cell.
    for (std::size_t k=0; k<nmat; ++k) {
      auto gkd = getDeformGradDev(nmat, k, solidx, ul);
      for (std::size_t i=0; i<3; ++i)
        for (std::size_t j=0; j<3; ++j) gk[i][j] = gkd[i][j];
      ur[energyIdx(nmat, k)] =
      mat_blk[k].template compute< EOS::totalenergy >(
        ur[densityIdx(nmat, k)], v1l, v2l, v3l,
        ul[volfracIdx(nmat, k)]*bp.fp, ul[volfracIdx(nmat, k)], gk );

      // material pressures
      ur[ncomp+pressureIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * bp.fp;
    }
  }
  // Otherwise, for supersonic outflow, all the characteristics are from
  // internal cell. Therefore, we calculate the ghost cell state using the
  // conservative variables from internal cell (which is what ur is
  // initialized to).
}

//! Back-pressure boundary state: deck pressure, internal-cell temperature
//! \details Mirror of inciter::back_pressure in BCFunctions.hpp. This is the
//!   only boundary state function that calls compute< EOS::temperature >.
//! \note As with inlet, the host lets the deformation gradient default to an
//!   empty std::array in both the temperature and totalenergy calls; a zero
//!   matrix is passed here to match. See the note on inletDev.
template< class MatBlkT, class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
backPressureDev( std::size_t ncomp,
                 std::size_t nmat,
                 std::size_t nstate,
                 const SolidxT& /*solidx*/,
                 const FnT& /*fn*/,
                 const StateT& ul,
                 StateT& ur,
                 const MatBlkT& mat_blk,
                 const BCParamsDev& bp )
{
  for (std::size_t i=0; i<nstate; ++i) ur[i] = ul[i];

  // Internal cell velocity components
  tk::real vl[3];
  for (std::size_t i=0; i<3; ++i)
    vl[i] = ul[ncomp+velocityIdx(nmat, i)];

  // see the note above: stands in for the host's defaulted empty defgrad
  const tk::real zerog[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};

  // The ghost cell state is calculated using the back pressure and other
  // quantities from the internal cell
  auto rhor_b(0.0);
  for (std::size_t k=0; k<nmat; ++k) {
    auto T_k = mat_blk[k].template compute< EOS::temperature >(
      ul[densityIdx(nmat, k)], vl[0], vl[1], vl[2], ul[energyIdx(nmat, k)],
      ul[volfracIdx(nmat, k)], zerog );
    auto rhor_k = mat_blk[k].template compute< EOS::density >( bp.fbp, T_k );
    ur[densityIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * rhor_k;
    ur[energyIdx(nmat, k)] =
    mat_blk[k].template compute< EOS::totalenergy >(
      ul[volfracIdx(nmat, k)]*rhor_k, vl[0], vl[1], vl[2],
      ul[volfracIdx(nmat, k)]*bp.fbp, ul[volfracIdx(nmat, k)], zerog );
    rhor_b += ur[densityIdx(nmat, k)];

    // material pressures
    ur[ncomp+pressureIdx(nmat, k)] = ul[volfracIdx(nmat, k)] * bp.fbp;
  }

  // bulk momentum
  for (std::size_t i=0; i<3; ++i)
    ur[momentumIdx(nmat,i)] = rhor_b * vl[i];
}

//! Dispatch to the boundary state function selected by kind
template< class MatBlkT, class FnT, class SolidxT, class StateT >
KOKKOS_INLINE_FUNCTION void
bcStateDev( BCKind kind,
            std::size_t ncomp,
            std::size_t nmat,
            std::size_t nstate,
            const SolidxT& solidx,
            const FnT& fn,
            const StateT& ul,
            StateT& ur,
            const MatBlkT& mat_blk,
            const BCParamsDev& bp )
{
  if (kind == BCKind::Symmetry)
    symmetryDev( ncomp, nmat, nstate, solidx, fn, ul, ur );
  else if (kind == BCKind::NoSlipWall)
    noslipwallDev( ncomp, nmat, nstate, solidx, fn, ul, ur );
  else if (kind == BCKind::Inlet)
    inletDev( ncomp, nmat, nstate, solidx, fn, ul, ur, mat_blk, bp );
  else if (kind == BCKind::Farfield)
    farfieldDev( ncomp, nmat, nstate, solidx, fn, ul, ur, mat_blk, bp );
  else if (kind == BCKind::BackPressure)
    backPressureDev( ncomp, nmat, nstate, solidx, fn, ul, ur, mat_blk, bp );
  else
    extrapolateDev( nstate, ul, ur );
}

} // inciter::

#endif // BCFunctionsDev_h
