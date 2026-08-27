// *****************************************************************************
/*!
  \file      src/PDE/Integrate/MultiMatTerms.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals of multi-material terms
     using DG methods
  \details   This file contains functionality for computing volume integrals of
     non-conservative and pressure relaxation terms that appear in the
     multi-material hydrodynamic equations, using the discontinuous Galerkin
     method for various orders of numerical representation.
*/
// *****************************************************************************
#ifndef MultiMatTerms_h
#define MultiMatTerms_h

#include "Basis.hpp"
#include "Fields.hpp"
#include "UnsMesh.hpp"
#include "EoS/EOS.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "Kokkos_Core.hpp"
#include "KokkosDevice.hpp"

using execution_space = Kokkos::DefaultExecutionSpace;
using memory_space = Kokkos::DefaultExecutionSpace::memory_space;

namespace tk {

using ncomp_t = tk::ncomp_t;

//! Compute volume integrals of non-conservative terms for multi-material DG
void
nonConservativeInt( const bool pref,
                    std::size_t nmat,
                    const std::vector< inciter::EOS >& mat_blk,
                    const std::size_t ndof,
                    const std::size_t rdof,
                    const std::size_t nelem,
                    const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord,
                    const Fields& geoElem,
                    const Fields& U,
                    const Fields& P,
                    const std::vector< std::vector< tk::real > >& riemannDeriv,
                    const std::vector< std::size_t >& ndofel,
                    Fields& R,
                    int intsharp );

// Persistent device-resident scratch buffers for nonConservativeInt_constP
// See KokkosDevice.hpp -> tk::KokkosDeviceViews for more detail
using nonConsvIntDeviceViews = KokkosDeviceViews;

// Kokkos implementation of nonConsvInt_constP
void
nonConservativeInt_constP( std::size_t nmat,
                           const std::vector< inciter::EOS >& mat_blk,
                           const std::size_t ndof,
                           const std::size_t rdof,
                           const std::size_t nelem,
                           const std::vector< std::size_t >& inpoel,
                           const UnsMesh::Coords& coord,
                           const Fields& geoElem,
                           const Fields& U,
                           const Fields& P,
                           const std::vector< std::vector< tk::real > >& riemannDeriv,
                           Fields& R,
                           int intsharp, 
                           nonConsvIntDeviceViews* dev = nullptr, bool prestaged = false );

// Kokkos device version of updateRhsNonCons
KOKKOS_INLINE_FUNCTION
void
updateRhsNonCons_device( ncomp_t ncomp,
                         const std::size_t nmat,
                         const std::size_t ndof,
                         const std::size_t ndof_el,
                         const tk::real wt,
                         const std::size_t r_nprop,
                         const std::size_t e,
                         const Kokkos::Array<tk::real, NDOF_MAX>& B,
                         const Kokkos::Array<Kokkos::Array<tk::real,NDOF_MAX>,NCOMP_MAX>& ncf,
                         Kokkos::View<real*, memory_space> R );

//! Update the rhs by adding the non-conservative term integrals
void
updateRhsNonCons( ncomp_t ncomp,
                const std::size_t nmat,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t e,
                const std::vector<tk::real>& B,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< std::vector< tk::real > >& ncf,
                Fields& R );

//! Compute volume integrals of non-conservative terms for multi-material FV
std::vector< tk::real >
nonConservativeIntFV(
  std::size_t nmat,
  const std::size_t rdof,
  const std::size_t e,
  const std::array< tk::real, 3 >& fn,
  const Fields& U,
  const Fields& P,
  const std::vector< tk::real >& var_riemann );

//! Compute volume integrals of pressure relaxation terms in multi-material DG
void
pressureRelaxationInt( const bool pref,
                       std::size_t nmat,
                       const std::vector< inciter::EOS >& mat_blk,
                       const std::size_t ndof,
                       const std::size_t rdof,
                       const std::size_t nelem,
                       const std::vector< std::size_t >& inpoel,
                       const UnsMesh::Coords& coord,
                       const Fields& geoElem,
                       const Fields& U,
                       const Fields& P,
                       const std::vector< std::size_t >& ndofel,
                       const tk::real ct,
                       Fields& R,
                       int intsharp );

//! Compute pressure relaxation integrals for const-order DG (not p-adaptive)
void
pressureRelaxationInt_constP( std::size_t nmat,
                              const std::vector< inciter::EOS >& mat_blk,
                              //const std::vector< tk::EOSDevice >& eosd,
                              const std::size_t ndof,
                              const std::size_t rdof,
                              const std::size_t nelem,
                              const std::vector< std::size_t >& inpoel,
                              const UnsMesh::Coords& coord,
                              const Fields& geoElem,
                              const Fields& U,
                              const Fields& P,
                              const tk::real ct,
                              Fields& R,
                              int intsharp,
                              KokkosDeviceViews* dev = nullptr,
                              bool prestaged = false );

//! Update the rhs by adding the pressure relaxation integrals
void
updateRhsPre(
  ncomp_t ncomp,
  const std::size_t ndof,
  const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector< tk::real >& B,
  std::vector< tk::real >& ncf,
  Fields& R );

//! Compute volume integrals of pressure relaxation terms in multi-material FV
void
pressureRelaxationIntFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  const tk::real ct,
  Fields& R );

//! Solve the reconstruct velocity used for volume fraction equation
std::vector< std::vector< tk::real > >
solvevriem( std::size_t nelem,
            const std::vector< std::vector< tk::real > >& vriem,
            const std::vector< std::vector< tk::real > >& riemannLoc );

//! Compute the riemann velociry at the interface
void evaluRiemann( ncomp_t ncomp,
                   const int e_left,
                   const int e_right,
                   const std::size_t nmat,
                   const std::vector< tk::real >& fl,
                   const std::array< tk::real, 3 >& fn,
                   const std::array< tk::real, 3 >& gp,
                   const std::array< std::vector< tk::real >, 2 >& state,
                   std::vector< std::vector< tk::real > >& vriem,
                   std::vector< std::vector< tk::real > >& riemannLoc );

//! Compute the flux-function for the multimaterial PDEs
std::vector< std::array< tk::real, 3 > >
fluxTerms(
  std::size_t ncomp,
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::vector< tk::real >& ugp );

//! Kokkos version of fluxTerms for experimental purposes
//! to avoid using function pointers
//! Streams one component's flux triple at a time into sink rather than filling
//! a large [NCOMP_MAX][3] array for the caller to walk afterwards
//! sink(c, f0, f1, f2) must apply the component-c contribution to R
//! Note FluxSink is a template param and NOT a std::function
template< typename FluxSink >
KOKKOS_INLINE_FUNCTION
void fluxTerms_multimat_kokkos(
  std::size_t ncomp,
  std::size_t nmat,
  Kokkos::View<const size_t*, memory_space>  solidx,
  //const std::vector< inciter::EOS >& /*mat_blk*/,
  Kokkos::Array<tk::real, NSTATE_MAX>& ugp, // this is state essentially
  Kokkos::Array<Kokkos::Array<Kokkos::Array<tk::real, 3>, 3>, NMAT_MAX>& g, 
  Kokkos::Array<Kokkos::Array<Kokkos::Array<tk::real, 3>, 3>, NMAT_MAX>& asig,
  Kokkos::Array<tk::real, NMAT_MAX>& al,
  //Kokkos::Array<Kokkos::Array<tk::real, 3>, NCOMP_MAX> & fl, 
  FluxSink&& sink,
  Kokkos::Array<tk::real, NMAT_MAX> & apk)
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::pressureIdx;
  using inciter::deformIdx;

  tk::real rho(0.0);
  for (std::size_t k=0; k<nmat; ++k)
    rho += ugp[densityIdx(nmat, k)];

  auto u = ugp[ncomp+velocityIdx(nmat,0)];
  auto v = ugp[ncomp+velocityIdx(nmat,1)];
  auto w = ugp[ncomp+velocityIdx(nmat,2)];
  //printf("u,v,w = %e, %e, %e\n", u, v, w);
 
  //? Have solid function is easy to do, its in PDE/MultiMat/MiscMultiMatFns.hpp //?DONE
  if (inciter::haveSolid(nmat, solidx))
  {
    Kokkos::Array<Kokkos::Array<real, 3>, 3> sig = {};
    for (std::size_t k=0; k<nmat; ++k)
    {
      al[k] = ugp[volfracIdx(nmat, k)];
      // inv deformation gradient and Cauchy stress tensors
      //? getDeformGrad, getCauchyStress seems fine, just need to pass in solidx to avoid using ginputDeck
      //? Try using subview, but doesnt seem to work? So, just pass g, asig as reference views :(
      inciter::getDeformGrad(nmat, k, solidx, ugp, g);
      inciter::getCauchyStress(nmat, k, ncomp, solidx, ugp, asig);
      for (std::size_t i=0; i<3; ++i) asig[k][i][i] -= ugp[ncomp+pressureIdx(nmat,k)];
      for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
          sig[i][j] += asig[k][i][j];
    }

    // conservative part of momentum flux
    for (std::size_t idir=0; idir<3; ++idir)
      sink( momentumIdx(nmat,idir),
            ugp[momentumIdx(nmat, idir)] * u - sig[0][idir],
            ugp[momentumIdx(nmat, idir)] * v - sig[1][idir],
            ugp[momentumIdx(nmat, idir)] * w - sig[2][idir] );

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      // skipped because it's zero

      // conservative part of material continuity flux
      sink( densityIdx(nmat,k), u*ugp[densityIdx(nmat, k)], 
                                v*ugp[densityIdx(nmat, k)],
                                w*ugp[densityIdx(nmat, k)] );

      // conservative part of material total-energy flux
      sink( energyIdx(nmat, k),
            u * ugp[energyIdx(nmat, k)]
              - u * asig[k][0][0] - v * asig[k][1][0] - w * asig[k][2][0],
            v * ugp[energyIdx(nmat, k)]
              - u * asig[k][0][1] - v * asig[k][1][1] - w * asig[k][2][1],
            w * ugp[energyIdx(nmat, k)]
              - u * asig[k][0][2] - v * asig[k][1][2] - w * asig[k][2][2] );

      // conservative part of material inverse deformation gradient
      // g_ij: \partial (g_il u_l) / \partial (x_j)
      if (solidx(k) > 0)
      {
        for (std::size_t i=0; i<3; ++i)
        {
          // No need to loop over j index so we can hoist this one
          const auto gu = u*g[k][i][0] + v*g[k][i][1] + w*g[k][i][2];
          for (std::size_t j=0; j<3; ++j)
          {
            sink( deformIdx(nmat, solidx(k), i, j), j==0?gu:0.0, j=1?gu:0.0, j=2?gu:0.0 );
          }
        }
      }
    }
  }
  else
  {
    tk::real p(0.0);
    for (std::size_t k=0; k<nmat; ++k)
    {
      apk[k] = ugp[ncomp+pressureIdx(nmat,k)];
      p += apk[k];
    }

    // conservative part of momentum flux
    // pressure contribution only on the diagonal (idir==direction)
    for (std::size_t idir=0; idir<3; ++idir)
      sink( momentumIdx(nmat, idir),
            ugp[momentumIdx(nmat, idir)] * u + (idir==0 ? p : 0.0),
            ugp[momentumIdx(nmat, idir)] * v + (idir==1 ? p : 0.0),
            ugp[momentumIdx(nmat, idir)] * w + (idir==2 ? p : 0.0) );

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      // skipped because it's zero

      // conservative part of material continuity flux
      sink( densityIdx(nmat, k), u * ugp[densityIdx(nmat, k)],
                                 v * ugp[densityIdx(nmat, k)],
                                 w * ugp[densityIdx(nmat, k)] );

      // conservative part of material total-energy flux
      auto hmat = ugp[energyIdx(nmat, k)] + apk[k];
      sink( energyIdx(nmat, k), u*hmat, v*hmat, w*hmat );
    }
  }
}
} // tk::

#endif // MultiMatTerms_h
