// *****************************************************************************
/*!
  \file      src/PDE/Integrate/MultiMatTerms.cpp
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

#include "QuinoaConfig.hpp"

#include "MultiMatTerms.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

// Lapacke forward declarations
extern "C" {

using lapack_int = long;

#define LAPACK_ROW_MAJOR 101

lapack_int LAPACKE_dsysv( int, char, lapack_int, lapack_int, double*,
    lapack_int, lapack_int*, double*, lapack_int );

}

// Local Kokkos helpers private to nonConsvIntConstP
namespace {
    using range_policy = Kokkos::RangePolicy<execution_space>;
    using UnManagedMem = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

    // Wrap host raw pointer in an unmanaged host view
    // Ensure all kernel entry points share one definition
    template<typename T>
    auto changeToView(T* object, std::size_t n) {
        Kokkos::View<T*, Kokkos::LayoutLeft, Kokkos::HostSpace, UnManagedMem>
          object_view(object, n);
        return object_view;
    }
}

namespace tk {

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
                    int intsharp )
// *****************************************************************************
//  Compute volume integrals for multi-material DG
//! \details This is called for multi-material DG, computing volume integrals of
//!   terms in the volume fraction and energy equations, which do not exist in
//!   the single-material flow formulation (for `CompFlow` DG). For further
//!   details see Pelanti, M., & Shyue, K. M. (2019). A numerical model for
//!   multiphase liquid–vapor–gas flows with interfaces and cavitation.
//!   International Journal of Multiphase Flow, 113, 208-230.
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface reconstruction indicator
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;
  using inciter::newSolidsAccFn;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > state(ncomp+nprim);

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    auto ng = tk::NGvol(ndofel[e]);

    // arrays for quadrature points
    std::array< std::vector< real >, 3 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    GaussQuadratureTet( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
            inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    // If an rDG method is set up (P0P1), then, currently we compute the P1
    // basis functions and solutions by default. This implies that P0P1 is
    // unsupported in the p-adaptive DG (PDG).
    std::size_t dof_el;
    if (rdof > ndof)
    {
      dof_el = rdof;
    }
    else
    {
      dof_el = ndofel[e];
    }

    // For multi-material p-adaptive simulation, when dofel = 1, p0p1 is
    // applied and ndof for solution evaluation should be 4
    if(dof_el == 1 && pref)
      dof_el = 4;

    // Compute the derivatives of basis function for second order terms
    std::array< std::vector<tk::real>, 3 > dBdx;
    for (std::size_t i=0; i<3; ++i) dBdx[i].resize( ndofel[e], 0 );
    if (ndofel[e] > 1)
      eval_dBdx_p1( ndofel[e], jacInv, dBdx );

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (ndofel[e] > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp], B );

      auto wt = wgp[igp] * geoElem(e, 0);

      evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
        rdof, nmat, e, dof_el, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);

      // get bulk properties
      tk::real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
          rhob += state[densityIdx(nmat, k)];

      // get the velocity vector
      std::array< tk::real, 3 > vel{{ state[ncomp+velocityIdx(nmat, 0)],
                                      state[ncomp+velocityIdx(nmat, 1)],
                                      state[ncomp+velocityIdx(nmat, 2)] }};

      std::vector< tk::real > ymat(nmat, 0.0);
      std::array< tk::real, 3 > dap{{0.0, 0.0, 0.0}};
      for (std::size_t k=0; k<nmat; ++k)
      {
        ymat[k] = state[densityIdx(nmat, k)]/rhob;

        std::size_t mark(3*k);
        if (solidx[k] > 0) mark = 3*nmat+ndof+3*(solidx[k]-1);

        for (std::size_t idir=0; idir<3; ++idir)
          dap[idir] += riemannDeriv[mark+idir][e];
      }

      // compute non-conservative terms
      std::vector< std::vector< tk::real > > ncf
        (ncomp, std::vector<tk::real>(ndofel[e],0.0));

      for (std::size_t idir=0; idir<3; ++idir)
        for(std::size_t idof=0; idof<ndofel[e]; ++idof)
          ncf[momentumIdx(nmat, idir)][idof] = 0.0;

      for (std::size_t k=0; k<nmat; ++k)
      {
        // evaluate non-conservative term for energy equation
        std::size_t mark(3*k);
        if (solidx[k] > 0) mark = 3*nmat+ndof+3*(solidx[k]-1);

        for(std::size_t idof=0; idof<ndofel[e]; ++idof)
        {
          ncf[densityIdx(nmat, k)][idof] = 0.0;

          for (std::size_t idir=0; idir<3; ++idir)
            ncf[energyIdx(nmat, k)][idof] -= vel[idir] * ( ymat[k]*dap[idir]
                                                  - riemannDeriv[mark+idir][e] );
        }

        // evaluate non-conservative terms for g equation
        if (solidx[k] > 0) {
          std::size_t nsld = inciter::numSolids(nmat, solidx);
          for (std::size_t idof=0; idof<ndofel[e]; ++idof)
            for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                for (std::size_t l=0; l<3; ++l)
                {
                  mark = 3*nmat+ndof+3*nsld+inciter::newSolidsAccFn(k,i,j,l);
                  ncf[deformIdx(nmat, solidx[k], i, j)][idof] -=
                    vel[l] * riemannDeriv[mark][e];
                }
        }

        // Evaluate non-conservative term for volume fraction equation:
        // Here we make an assumption that the derivative of Riemann velocity
        // times the basis function is constant. Therefore, when P0P1/DGP1/DGP2
        // are used for constant velocity problems, the discretization is
        // consistent. However, for a general problem with varying velocity,
        // there will be errors since the said derivative is not constant.
        // A discretization that solves this issue has not been implemented yet.
        // Nevertheless, this does not affect high-order accuracy in
        // single material regions for problems with sharp interfaces. Since
        // volume fractions are nearly constant in such regions, using
        // high-order for volume fractions does not show any benefits over
        // THINC. Therefore, for such problems, we only use FV for the volume
        // fractions, and these non-conservative high-order terms do not need
        // to be computed.
        // In summary, high-order discretization for non-conservative terms in
        // volume fraction equations is avoided for sharp interface problems.
        if (ndof <= 4 || intsharp == 1) {
          for(std::size_t idof=0; idof<ndofel[e]; ++idof)
            ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
                                           * riemannDeriv[3*nmat+idof][e];
        } else if (intsharp == 0) {     // If DGP2 without THINC
          // DGP2 is discretized differently than DGP1/FV to guarantee 3rd order
          // convergence for the testcases with uniform and constant velocity.

          // P0 contributions for all equations
          for(std::size_t idof=0; idof<ndof; ++idof)
          ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
                                         * riemannDeriv[3*nmat][e] * B[idof];
          // High order contributions
          for(std::size_t idof=1; idof<ndof; ++idof)
            for(std::size_t idir=0; idir<3; ++idir)
            ncf[volfracIdx(nmat, k)][idof] += state[volfracIdx(nmat, k)]
                                            * vel[idir] * dBdx[idir][idof];
        }
      }

      updateRhsNonCons( ncomp, nmat, ndof, ndofel[e], wt, e, B, dBdx, ncf, R );
    }
  }
}

void
nonConservativeInt_constP(
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
  const std::vector<std::vector<tk::real>>& riemannDeriv,
  Fields& R,
  int intsharp,
  nonConsvIntDeviceViews* dev, bool prestaged ) //added this
// *****************************************************************************
//  Compute volume integrals for multi-material DG (const-order, not p-adaptive)
//! \details This is called for multi-material DG, computing volume integrals of
//!   terms in the volume fraction and energy equations, which do not exist in
//!   the single-material flow formulation (for `CompFlow` DG). For further
//!   details see Pelanti, M., & Shyue, K. M. (2019).
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface reconstruction indicator
// *****************************************************************************
{
  Kokkos::Profiling::pushRegion("nonConsvInt");

  using inciter::volfracIdx;
  using inciter::densityIdx;
  //using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;
  //using inciter::newSolidsAccFn;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  // Interface-compression parameter read on host
  // g_inputdeck is not addressable from the device
  auto bparam = inciter::g_inputdeck.get<tag::multimat,tag::intsharp_param>();

  const std::size_t ncomp = U.nprop()/rdof;
  const std::size_t nprim = P.nprop()/rdof;

  const std::size_t m_nprop = U.nprop();
  const std::size_t p_nprop = P.nprop();
  const std::size_t r_nprop = R.nprop();
  const std::size_t geo_nprop = geoElem.nprop();

  // Fail loudly on the hoest if the config overruns fixed-size device scratch arrays
  // Otherwise we are corrupting thread-private memory
  checkKokkosCaps( nmat, ndof, rdof, ncomp, nprim );

  // Need to evaluate numSolids() here on host before capturing it into the kernel
  const std::size_t nsld = inciter::numSolids(nmat, solidx);

  // Quadrature points
  // We only need to compute this once for constant P
  const auto ng = tk::NGvol(ndof);

  // Persistent device buffers
  // Prevents redundant cudaMallocs and cudaFrees per call
  nonConsvIntDeviceViews local_dev;
  nonConsvIntDeviceViews& dv = dev ? *dev : local_dev;

  // Variables go here
  // First check whether device already holds this partition's mesh data at current meshgen state
  // Refer to tk::meshResident() for more details
  const bool mesh_ok = meshResident( dv, inpoel, coord, geoElem, nelem, nmat );

  // Solidx
  auto solidx_h_view = changeToView(solidx.data(), nmat);
  if (ensureDeviceCapacity(dv.solidx, "nconsv_solidx_d_view", nmat) || !mesh_ok)
    Kokkos::deep_copy(dv.solidx, solidx_h_view);

  // Inpoel
  const std::size_t inpoel_size = inpoel.size();
  auto inpoel_h_view = changeToView(inpoel.data(), inpoel_size);
  if (ensureDeviceCapacity(dv.inpoel, "nconsv_inpoel_d_view", inpoel_size) || !mesh_ok)
    Kokkos::deep_copy(dv.inpoel, inpoel_h_view);

  // Transfer coord (nodal coordinates)
  size_t coordx_size = coord[0].size();
  auto cx_h_view = changeToView(coord[0].data(), coordx_size);
  if (ensureDeviceCapacity(dv.cx, "nconsv_cx_d_view", coordx_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cx, cx_h_view);

  size_t coordy_size = coord[1].size();
  auto cy_h_view = changeToView(coord[1].data(), coordy_size);
  if (ensureDeviceCapacity(dv.cy, "nconsv_cy_d_view", coordy_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cy, cy_h_view);
  
  size_t coordz_size = coord[2].size();
  auto cz_h_view = changeToView(coord[2].data(), coordz_size);
  if (ensureDeviceCapacity(dv.cz, "nconsv_cz_d_view", coordz_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cz, cz_h_view);

  // geoElem, U, P transfer
  size_t geoElem_size = geoElem.getSize();
  auto geoElem_h_view = changeToView(geoElem.getPointer(), geoElem_size);
  if (ensureDeviceCapacity(dv.geoElem, "nconsv_geoElem_d_view", geoElem_size) || !mesh_ok)
    Kokkos::deep_copy(dv.geoElem, geoElem_h_view);

  /*
  // U,P change per call so always upload
  size_t P_size = P.getSize();
  ensureDeviceCapacity(dv.P, "nconsv_P_d_view", P_size);
  auto P_h_view = changeToView(P.getPointer(), P_size);
  Kokkos::deep_copy(dv.P, P_h_view);

  size_t U_size = U.getSize();
  ensureDeviceCapacity(dv.U, "nconsv_U_d_view", U_size);
  auto U_h_view = changeToView(U.getPointer(), U_size);
  Kokkos::deep_copy(dv.U, U_h_view);
  */

  // Up, P, R change per call so always upload
  // Staged thru pagelock memory and queued on default exec space instance (async copies)
  // Kernel runs on the same instance (stream) so its ordered after without needing to fence
  auto exec = Kokkos::DefaultExecutionSpace();

  const std::size_t R_size = R.getSize();
  if(!prestaged){
    uploadStaged( exec, dv.R, dv.stage_R, R.getPointer(), R_size, "nconsv_R_d_view" );

    const std::size_t P_size = P.getSize();
    uploadStaged( exec, dv.P, dv.stage_P, P.getPointer(), P_size, "nconsv_P_d_view" );

    const std::size_t U_size = U.getSize();
    uploadStaged( exec, dv.U, dv.stage_U, U.getPointer(), U_size, "nconsv_U_d_view" );
  }

  // RiemannDeriv only consumed here, nothing on the host runs between surf ints and this kernel
  // So, there's nothing to hide it behind and we want to stage it locally
  const std::size_t rd_nrow = riemannDeriv.size();
  const std::size_t rd_ncol = rd_nrow ? riemannDeriv[0].size() : 0;
  const std::size_t rd_size = rd_nrow*rd_ncol;
  std::vector<tk::real> rd_flat(rd_size);
  for (std::size_t row=0; row<rd_nrow; ++row){
    for (std::size_t col=0; col<rd_ncol; ++col){
      rd_flat[row*rd_ncol + col] = riemannDeriv[row][col];
    }
  }
  uploadStaged( exec, dv.riemannDeriv, dv.stage_rd, rd_flat.data(), rd_size, "nconsv_riemannDeriv_d_view" );

  // Shallow copies of view handle
  // Does not touch device memory, just gives kernel below local names and captures plain views by value into KOKKOS_LAMBDA
  // If volIntDeviceViews struct (or a ref to it) is captured it will be invalid on device
  auto solidx_d_view = dv.solidx;
  auto inpoel_d_view = dv.inpoel;
  auto cx_d_view = dv.cx;
  auto cy_d_view = dv.cy;
  auto cz_d_view = dv.cz;
  auto geoElem_d_view = dv.geoElem;
  auto P_d_view = dv.P;
  auto U_d_view = dv.U;
  auto rd_d_view = dv.riemannDeriv;
  auto R_d_view = dv.R;

  // Quadrature points can be hoisted out
  // Because P is constant
  Kokkos::Array<Kokkos::Array<real,NQUAD_MAX>, 3> coordgp = {};
  Kokkos::Array<real,NQUAD_MAX> wgp = {};
  GaussQuadratureTet(ng,coordgp,wgp);

  // compute volume integrals
  Kokkos::parallel_for("nonConservativeInt_kernel",range_policy(0,nelem),KOKKOS_LAMBDA(const std::size_t e)
  {
    // Extract elem coords
    Kokkos::Array<Kokkos::Array<real,3>,4> coordel;
    for (int i=0; i<4; ++i){
      coordel[i][0] = cx_d_view(inpoel_d_view(4*e+i));
      coordel[i][1] = cy_d_view(inpoel_d_view(4*e+i));
      coordel[i][2] = cz_d_view(inpoel_d_view(4*e+i));
    }

    auto jacInv =
            inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    // Compute the derivatives of basis function for second order terms
    // eval_dBdx_p1( ndof, jacInv, dBdx );
    Kokkos::Array<real, NDOF_MAX> B = {};
    Kokkos::Array<Kokkos::Array<real, NDOF_MAX>, 3> dBdx = {};
    eval_dBdx_p1(ndof, jacInv, dBdx);

    Kokkos::Array<real, NSTATE_MAX> state = {};

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (ndof > 4) eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the basis function
      eval_basis( rdof, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp], B );

      //auto wt = wgp[igp] * geoElem(e, 0);
      auto wt = wgp[igp] * geoElem_d_view(e*geo_nprop);

      evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
        //rdof, nmat, e, rdof, inpoel, coord, geoElem,
        //{{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);
        rdof,nmat,e,rdof,m_nprop,p_nprop,geo_nprop,
        bparam,solidx_d_view,inpoel_d_view,
        cx_d_view, cy_d_view, cz_d_view, geoElem_d_view,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U_d_view, P_d_view, state);

      // get bulk properties
      // tk::real rhob(0.0);
      real rhob(0.0);
      for (std::size_t k=0; k<nmat; ++k)
          rhob += state[densityIdx(nmat, k)];

      // get the velocity vector
      //std::array< tk::real, 3 > vel{{ state[ncomp+velocityIdx(nmat, 0)],
      Kokkos::Array< real, 3 > vel{{ state[ncomp+velocityIdx(nmat, 0)],
                                      state[ncomp+velocityIdx(nmat, 1)],
                                      state[ncomp+velocityIdx(nmat, 2)] }};

      //std::vector< tk::real > ymat(nmat, 0.0);
      //std::array< tk::real, 3 > dap{{0.0, 0.0, 0.0}};
      Kokkos::Array<real, NMAT_MAX> ymat = {};
      Kokkos::Array<real,3> dap{{0.0,0.0,0.0}};
      for (std::size_t k=0; k<nmat; ++k)
      {
        ymat[k] = state[densityIdx(nmat, k)]/rhob;

        std::size_t mark(3*k);
        //if (solidx[k] > 0) mark = 3*nmat+ndof+3*(solidx[k]-1);
        if(solidx_d_view(k)>0) mark = 3*nmat+ndof+3*(solidx_d_view(k)-1);

        for (std::size_t idir=0; idir<3; ++idir)
          //dap[idir] += riemannDeriv[mark+idir][e];
          dap[idir] += rd_d_view((mark+idir)*rd_ncol + e);
      }

      /*
      // compute non-conservative terms
      std::vector< std::vector< tk::real > > ncf
        (ncomp, std::vector<tk::real>(ndof,0.0));

      for (std::size_t idir=0; idir<3; ++idir)
        for(std::size_t idof=0; idof<ndof; ++idof)
          ncf[momentumIdx(nmat, idir)][idof] = 0.0;
      */
      // zero-initialize and compute nonconsv term
      Kokkos::Array<Kokkos::Array<real, NDOF_MAX>, NCOMP_MAX> ncf = {};

      for (std::size_t k=0; k<nmat; ++k)
      {
        // evaluate non-conservative term for energy equation
        std::size_t mark(3*k);
        if (solidx_d_view(k) > 0) mark = 3*nmat+ndof+3*(solidx_d_view(k)-1);
        //if (solidx[k] > 0) mark = 3*nmat+ndof+3*(solidx[k]-1);

        for(std::size_t idof=0; idof<ndof; ++idof)
        {
          ncf[densityIdx(nmat, k)][idof] = 0.0;

          for (std::size_t idir=0; idir<3; ++idir)
            ncf[energyIdx(nmat, k)][idof] -= vel[idir] * ( ymat[k]*dap[idir]
            //                                      - riemannDeriv[mark+idir][e] );
                                                  - rd_d_view((mark+idir)*rd_ncol + e) );
        }

        // evaluate non-conservative terms for g equation
        //if (solidx[k] > 0) {
        //  std::size_t nsld = inciter::numSolids(nmat, solidx);
        if (solidx_d_view(k) > 0) {
          for (std::size_t idof=0; idof<ndof; ++idof)
            for (std::size_t i=0; i<3; ++i)
              for (std::size_t j=0; j<3; ++j)
                for (std::size_t l=0; l<3; ++l)
                {
                  //mark = 3*nmat+ndof+3*nsld+inciter::newSolidsAccFn(k,i,j,l);
                  std::size_t gmark = 3*nmat+ndof+3*nsld + inciter::newSolidsAccFn(k,i,j,l);
                  ncf[deformIdx(nmat, solidx_d_view(k),i,j)][idof] -=
                    vel[l]*rd_d_view(gmark*rd_ncol+e);
                  //ncf[deformIdx(nmat, solidx[k], i, j)][idof] -=
                  //  vel[l] * riemannDeriv[mark][e];
                }
        }

        // Evaluate non-conservative term for volume fraction equation.
        // (See original notes on the constant-Riemann-velocity assumption and
        // the special DGP2 discretization.)
        if (ndof <= 4 || intsharp == 1) {
          for(std::size_t idof=0; idof<ndof; ++idof)
            ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
                                           * rd_d_view((3*nmat+idof)*rd_ncol+e);
            //                               * riemannDeriv[3*nmat+idof][e];
        } else if (intsharp == 0) {     // If DGP2 without THINC
          // DGP2 is discretized differently than DGP1/FV to guarantee 3rd order
          // convergence for testcases with uniform and constant velocity.

          // P0 contributions for all equations
          for(std::size_t idof=0; idof<ndof; ++idof)
            ncf[volfracIdx(nmat, k)][idof] = state[volfracIdx(nmat, k)]
          //                                 * riemannDeriv[3*nmat][e] * B[idof];
                                           * rd_d_view((3*nmat)*rd_ncol + e) * B[idof];
          // High order contributions
          for(std::size_t idof=1; idof<ndof; ++idof)
            for(std::size_t idir=0; idir<3; ++idir)
              ncf[volfracIdx(nmat, k)][idof] += state[volfracIdx(nmat, k)]
                                              * vel[idir] * dBdx[idir][idof];
        }
      }

      tk::updateRhsNonCons_device(ncomp,nmat,ndof,ndof,wt,r_nprop,e,B,ncf,R_d_view);
      //updateRhsNonCons( ncomp, nmat, ndof, ndof, wt, e, B, dBdx, ncf, R );
    }
  }); //end Kokkos::parallel_for

  //Kokkos::deep_copy(R_h_view, R_d_view);

  // Queued on same instance as kernel so ordered after it
  // Fence before copying out of pinned buffer
  if (!prestaged)
    downloadStaged( exec, R.getPointerNonConst(), dv.stage_R, R_d_view, R_size, "nconsv_R_d_view" ); 
 
  Kokkos::Profiling::popRegion();
}

//! Kokkos device version of updateRhsNonCons
KOKKOS_INLINE_FUNCTION
void updateRhsNonCons_device(
  ncomp_t ncomp,
  const std::size_t nmat,
  const std::size_t ndof,
  const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t r_nprop,
  const std::size_t e,
  //const Kokkos::Array<tk::real, 10>& B,
  //const Kokkos::Array<Kokkos::Array<tk::real, 10>, 50>& ncf,
  const Kokkos::Array<tk::real, NDOF_MAX>& B,
  const Kokkos::Array<Kokkos::Array<tk::real, NDOF_MAX>, NCOMP_MAX>& ncf,
  Kokkos::View<real*, memory_space> R)
// *****************************************************************************
//  Update the rhs by adding the non-conservative term integrals (Kokkos ver)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Basis function evaluated at local quadrature point
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] ncf Vector of non-conservative terms
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::energyIdx;
  using inciter::volfracDofIdx;
  using inciter::energyDofIdx;

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark=c*ndof;
    R(e*r_nprop+mark) += wt*ncf[c][0];
  }

  if(ndof_el>1)
  {
    // Update rhs with distributions from vol fraction and energy eqns
    for (std::size_t k=0; k<nmat; ++k)
    {
      for (std::size_t idof=1; idof<ndof_el; ++idof)
      {
        R(e * r_nprop + volfracDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[volfracIdx(nmat,k)][idof];
        R(e * r_nprop + energyDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[energyIdx(nmat,k)][idof] * B[idof];
        // Note: High order non-conservative g terms not implemented yet!
      }
    }
  }
}

void
updateRhsNonCons(
  ncomp_t ncomp,
  const std::size_t nmat,
  const std::size_t ndof,
  [[maybe_unused]] const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector<tk::real>& B,
  [[maybe_unused]] const std::array< std::vector<tk::real>, 3 >& dBdx,
  const std::vector< std::vector< tk::real > >& ncf,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the non-conservative term integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Basis function evaluated at local quadrature point
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] ncf Vector of non-conservative terms
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::energyIdx;
  using inciter::deformIdx;
  using inciter::volfracDofIdx;
  using inciter::energyDofIdx;
  using inciter::deformDofIdx;

  //Assert( dBdx[0].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[1].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[2].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( ncf.size() == ncomp,
  //        "Size mismatch for non-conservative term" );
  Assert( ncf.size() == ncomp, "Size mismatch for non-conservative term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark) += wt * ncf[c][0];
  }

  if( ndof_el > 1)
  {
    // Update rhs with distributions from volume fraction and energy equations
    for (std::size_t k=0; k<nmat; ++k)
    {
      for(std::size_t idof = 1; idof < ndof_el; idof++)
      {
        R(e, volfracDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[volfracIdx(nmat,k)][idof];
        R(e, energyDofIdx(nmat,k,ndof,idof)) +=
          wt * ncf[energyIdx(nmat,k)][idof] * B[idof];
        // Note: High order non-conservative g terms not implemented yet!
      }
    }
  }
}

std::vector< tk::real >
nonConservativeIntFV(
  std::size_t nmat,
  const std::size_t rdof,
  const std::size_t e,
  const std::array< tk::real, 3 >& fn,
  const Fields& U,
  const Fields& P,
  const std::vector< tk::real >& var_riemann )
// *****************************************************************************
//  Compute integrals of non-conservative terms for multi-material FV
//! \details This is called for multi-material FV, computing integrals of
//!   terms in the volume fraction and energy equations, which do not exist in
//!   the single-material flow formulation (for `CompFlow`). For further
//!   details see Pelanti, M., & Shyue, K. M. (2019). A numerical model for
//!   multiphase liquid–vapor–gas flows with interfaces and cavitation.
//!   International Journal of Multiphase Flow, 113, 208-230.
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] e Element for which contribution is to be calculated
//! \param[in] fn Face normal
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] var_riemann Riemann-values of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::volfracDofIdx;
  using inciter::densityDofIdx;
  using inciter::velocityDofIdx;

  auto ncomp = U.nprop()/rdof;

  // get bulk properties
  tk::real rhob(0.0), p_face(0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    rhob += U(e, densityDofIdx(nmat,k,rdof,0));
    p_face += var_riemann[k];
  }

  std::array< tk::real, 3 > vel{{ P(e, velocityDofIdx(nmat,0,rdof,0)),
                                  P(e, velocityDofIdx(nmat,1,rdof,0)),
                                  P(e, velocityDofIdx(nmat,2,rdof,0)) }};
  auto v_dot_n = tk::dot(vel, fn);

  // compute non-conservative terms
  auto v_riem = var_riemann[nmat];
  std::vector< tk::real > ncf(ncomp, 0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    auto ymat = U(e, densityDofIdx(nmat,k,rdof,0))/rhob;

    // evaluate non-conservative term for energy equation
    ncf[energyIdx(nmat, k)] -= v_dot_n * ( ymat*p_face - var_riemann[k] );

    // evaluate non-conservative term for volume fraction equation
    ncf[volfracIdx(nmat, k)] = U(e, volfracDofIdx(nmat,k,rdof,0))
      * v_riem;
  }

  return ncf;
}

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
                       int intsharp )
// *****************************************************************************
//  Compute volume integrals of pressure relaxation terms in multi-material DG
//! \details This is called for multi-material DG to compute volume integrals of
//!   finite pressure relaxation terms in the volume fraction and energy
//!   equations, which do not exist in the single-material flow formulation (for
//!   `CompFlow` DG). For further details see Dobrev, V. A., Kolev, T. V.,
//!   Rieben, R. N., & Tomov, V. Z. (2016). Multi‐material closure model for
//!   high‐order finite element Lagrangian hydrodynamics. International Journal
//!   for Numerical Methods in Fluids, 82(10), 689-706.
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in] ct Pressure relaxation time-scale for this system
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface reconstruction indicator
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::pressureIdx;
  using inciter::velocityIdx;
  using inciter::deformIdx;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > state(ncomp+nprim);

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    auto dx = geoElem(e,4)/2.0;
    auto ng = NGvol(ndofel[e]);

    // arrays for quadrature points
    std::array< std::vector< real >, 3 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    coordgp[2].resize( ng );
    wgp.resize( ng );

    GaussQuadratureTet( ng, coordgp, wgp );

    // Compute the derivatives of basis function for DG(P1)
    std::array< std::vector<real>, 3 > dBdx;

    // If an rDG method is set up (P0P1), then, currently we compute the P1
    // basis functions and solutions by default. This implies that P0P1 is
    // unsupported in the p-adaptive DG (PDG).
    std::size_t dof_el;
    if (rdof > ndof)
    {
      dof_el = rdof;
    }
    else
    {
      dof_el = ndofel[e];
    }

    // For multi-material p-adaptive simulation, when dofel = 1, p0p1 is applied
    // and ndof for solution evaluation should be 4
    if(dof_el == 1 && pref)
      dof_el = 4;

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp], B );

      auto wt = wgp[igp] * geoElem(e, 0);

      evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
        rdof, nmat, e, dof_el, inpoel, coord, geoElem,
        {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);

      // get bulk pressures and bulk modulii
      real pb(0.0), nume(0.0), deno(0.0), trelax(0.0);
      std::vector< real > apmat(nmat, 0.0), kmat(nmat, 0.0);
      std::vector< int > do_relax(nmat, 1);
      bool is_relax(false);
      for (std::size_t k=0; k<nmat; ++k)
      {
        real arhomat = state[densityIdx(nmat, k)];
        real alphamat = state[volfracIdx(nmat, k)];
        apmat[k] = state[ncomp+pressureIdx(nmat, k)];
        real amat = 0.0;
        bool include_solid(true);
        if (solidx[k] > 0 && apmat[k] < 1e3*alphamat) include_solid = false;
        if (include_solid && alphamat >= inciter::volfracPRelaxLim()) {
            amat = mat_blk[k].compute< inciter::EOS::soundspeed >( arhomat,
              apmat[k], alphamat, k );
          kmat[k] = arhomat * amat * amat;
          pb += apmat[k];

          // relaxation parameters
          trelax = std::max(trelax, ct*dx/amat);
          nume += alphamat * apmat[k] / kmat[k];
          deno += alphamat * alphamat / kmat[k];

          is_relax = true;
        }
        else do_relax[k] = 0;
      }
      real p_relax(0.0);
      if (is_relax) p_relax = nume/deno;

      // compute pressure relaxation terms
      std::vector< real > s_prelax(ncomp, 0.0);
      for (std::size_t k=0; k<nmat; ++k)
      {
        // only perform prelax on existing quantities
        if (do_relax[k] == 1) {
          auto s_alpha = (apmat[k]-p_relax*state[volfracIdx(nmat, k)])
            * (state[volfracIdx(nmat, k)]/kmat[k]) / trelax;
          s_prelax[volfracIdx(nmat, k)] = s_alpha;
          s_prelax[energyIdx(nmat, k)] = - pb*s_alpha;
        }
      }

      updateRhsPre( ncomp, ndof, ndofel[e], wt, e, B, s_prelax, R );
    }
  }
}

void
updateRhsPre(
  ncomp_t ncomp,
  const std::size_t ndof,
  [[maybe_unused]] const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector< tk::real >& B,
  std::vector< tk::real >& ncf,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the pressure relaxation integrals
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Basis function evaluated at local quadrature point
//! \param[in] ncf Vector of non-conservative terms
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  //Assert( dBdx[0].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[1].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( dBdx[2].size() == ndof_el,
  //        "Size mismatch for basis function derivatives" );
  //Assert( ncf.size() == ncomp,
  //        "Size mismatch for non-conservative term" );
  Assert( ncf.size() == ncomp, "Size mismatch for non-conservative term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    for(std::size_t idof = 0; idof < ndof; idof++)
      R(e, mark+idof) += wt * ncf[c] * B[idof];
  }
}

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
  Fields& R )
// *****************************************************************************
//  Compute volume integrals of pressure relaxation terms in multi-material FV
//! \details This is called for multi-material FV to compute volume integrals of
//!   finite pressure relaxation terms in the volume fraction and energy
//!   equations, which do not exist in the single-material flow formulation (for
//!   `CompFlow`). For further details see Dobrev, V. A., Kolev, T. V.,
//!   Rieben, R. N., & Tomov, V. Z. (2016). Multi‐material closure model for
//!   high‐order finite element Lagrangian hydrodynamics. International Journal
//!   for Numerical Methods in Fluids, 82(10), 689-706.
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Total number of elements
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitive quantities at recent time step
//! \param[in] ct Pressure relaxation time-scale for this system
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::volfracDofIdx;
  using inciter::energyDofIdx;
  using inciter::pressureIdx;
  using inciter::densityIdx;

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< real > apmat(nmat, 0.0), kmat(nmat, 0.0);
  std::vector< int > do_relax(nmat, 1);

  // Compute the basis function
  std::vector< tk::real > B(rdof, 0.0);
  std::vector< tk::real > state(ncomp+nprim);
  B[0] = 1.0;

  // compute volume integrals
  for (std::size_t e=0; e<nelem; ++e)
  {
    auto dx = geoElem(e,4)/2.0;

    evalPolynomialSol(mat_blk, 0, ncomp, nprim,
      rdof, nmat, e, rdof, inpoel, coord, geoElem,
      {{0.25, 0.25, 0.25}}, B, U, P, state);

    // get bulk pressures and bulk modulii
    real pb(0.0), nume(0.0), deno(0.0), trelax(0.0);
    bool is_relax(false);
    for (std::size_t k=0; k<nmat; ++k)
    {
      do_relax[k] = 1;
      real arhomat = state[densityIdx(nmat, k)];
      real alphamat = state[volfracIdx(nmat, k)];
      if (alphamat >= inciter::volfracPRelaxLim()) {
        apmat[k] = state[ncomp+pressureIdx(nmat, k)];
        real amat = mat_blk[k].compute< inciter::EOS::soundspeed >( arhomat,
          apmat[k], alphamat, k );
        kmat[k] = arhomat * amat * amat;
        pb += apmat[k];

        // relaxation parameters
        trelax = std::max(trelax, ct*dx/amat);
        nume += alphamat * apmat[k] / kmat[k];
        deno += alphamat * alphamat / kmat[k];

        is_relax = true;
      }
      else do_relax[k] = 0;
    }
    real p_relax(0.0);
    if (is_relax) p_relax = nume/deno;

    // compute pressure relaxation terms
    for (std::size_t k=0; k<nmat; ++k)
    {
      // only perform prelax on existing quantities
      if (do_relax[k] == 1) {
        auto s_alpha = (apmat[k]-p_relax*state[volfracIdx(nmat, k)])
          * (state[volfracIdx(nmat, k)]/kmat[k]) / trelax;

        R(e, volfracDofIdx(nmat,k,1,0)) += geoElem(e,0) * s_alpha;
        R(e, energyDofIdx(nmat,k,1,0)) += geoElem(e,0) * (-pb*s_alpha);
      }
    }
  }
}

std::vector< std::vector< tk::real > >
solvevriem( std::size_t nelem,
            const std::vector< std::vector< tk::real > >& vriem,
            const std::vector< std::vector< tk::real > >& riemannLoc )
// *****************************************************************************
//  Solve the reconstruct velocity used for volume fraction equation by
//  Least square method
//! \param[in] nelem Numer of elements
//! \param[in,out] vriem Vector of the riemann velocity
//! \param[in,out] riemannLoc Vector of coordinates where Riemann velocity data
//!   is available
//! \return Vector of Riemann velocity polynomial solution
// *****************************************************************************
{
  std::vector< std::vector< tk::real > >
    vriempoly( nelem, std::vector<tk::real>(12,0.0) );

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Use the normal method to construct the linear system A^T * A * x = u
    auto numgp = riemannLoc[e].size()/3;
    std::vector< std::vector< tk::real > > A(numgp,
                                             std::vector<tk::real>(4, 1.0));

    for(std::size_t k = 0; k < numgp; k++)
    {
      auto mark = k * 3;
      A[k][1] = riemannLoc[e][mark];
      A[k][2] = riemannLoc[e][mark+1];
      A[k][3] = riemannLoc[e][mark+2];
    }

    for(std::size_t idir = 0; idir < 3; idir++)
    {
      double AA_T[4*4], u[4];

      for(std::size_t i = 0; i < 4; i++)
        for(std::size_t j = 0; j < 4; j++)
        {
          auto id = 4 * i + j;
          AA_T[id] = 0;
          for(std::size_t k = 0; k < numgp; k++)
            AA_T[id] += A[k][i] * A[k][j];
        }

      std::vector<tk::real> vel(numgp, 1.0);
      for(std::size_t k = 0; k < numgp; k++)
      {
        auto mark = k * 3 + idir;
        vel[k] = vriem[e][mark];
      }
      for(std::size_t k = 0; k < 4; k++)
      {
        u[k] = 0;
        for(std::size_t i = 0; i < numgp; i++)
          u[k] += A[i][k] * vel[i];
      }
 
      lapack_int IPIV[4];
      #ifndef NDEBUG
      lapack_int info =
      #endif
        LAPACKE_dsysv( LAPACK_ROW_MAJOR, 'U', 4, 1, AA_T, 4, IPIV, u, 1 );
      Assert( info == 0, "Error in linear system solver" );

      auto idirmark = idir * 4;
      for(std::size_t k = 0; k < 4; k++)
        vriempoly[e][idirmark+k] = u[k];
    }
  }
  return vriempoly;
}

void evaluRiemann( ncomp_t ncomp,
                   const int e_left,
                   const int e_right,
                   const std::size_t nmat,
                   const std::vector< tk::real >& fl,
                   const std::array< tk::real, 3 >& fn,
                   const std::array< tk::real, 3 >& gp,
                   const std::array< std::vector< tk::real >, 2 >& state,
                   std::vector< std::vector< tk::real > >& vriem,
                   std::vector< std::vector< tk::real > >& riemannLoc )
// *****************************************************************************
//  Compute the riemann velocity at the interface
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] e_left Index for the left element
//! \param[in] e_right Index for the right element
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] fn Face/Surface normal
//! \param[in] gp Gauss points coordinates
//! \param[in] fl Surface flux
//! \param[in] state Vector of state variables for left and right side
//! \param[in,out] vriem Vector of the riemann velocity
//! \param[in,out] riemannLoc Vector of coordinates where Riemann velocity data
//!   is available
// *****************************************************************************
{
  using inciter::densityIdx;
  using inciter::momentumIdx;

  std::size_t el(0), er(0);
  el = static_cast< std::size_t >(e_left);
  if(e_right != -1)
    er = static_cast< std::size_t >(e_right);

  riemannLoc[el].push_back( gp[0] );
  riemannLoc[el].push_back( gp[1] );
  riemannLoc[el].push_back( gp[2] );

  if(e_right != -1)
  {
    riemannLoc[er].push_back( gp[0] );
    riemannLoc[er].push_back( gp[1] );
    riemannLoc[er].push_back( gp[2] );
  }

  tk::real rhobl(0.0), rhobr(0.0);
  for (std::size_t k=0; k<nmat; ++k)
  {
    rhobl += state[0][densityIdx(nmat, k)];
    rhobr += state[1][densityIdx(nmat, k)];
  }

  auto ul = state[0][momentumIdx(nmat, 0)] / rhobl;
  auto vl = state[0][momentumIdx(nmat, 1)] / rhobl;
  auto wl = state[0][momentumIdx(nmat, 2)] / rhobl;

  auto ur = state[1][momentumIdx(nmat, 0)] / rhobr;
  auto vr = state[1][momentumIdx(nmat, 1)] / rhobr;
  auto wr = state[1][momentumIdx(nmat, 2)] / rhobr;

  // Compute the normal velocities from left and right cells
  auto vnl = ul * fn[0] + vl * fn[1] + wl * fn[2];
  auto vnr = ur * fn[0] + vr * fn[1] + wr * fn[2];

  // The interface velocity is evaluated by adding the normal velocity which
  // is taken from the Riemann solver and the tangential velocity which is
  // evaluated as an average of the left and right cells
  auto urie = 0.5 * ((ul + ur) - fn[0] * (vnl + vnr)) + fl[ncomp+nmat] * fn[0];
  auto vrie = 0.5 * ((vl + vr) - fn[1] * (vnl + vnr)) + fl[ncomp+nmat] * fn[1];
  auto wrie = 0.5 * ((wl + wr) - fn[2] * (vnl + vnr)) + fl[ncomp+nmat] * fn[2];

  vriem[el].push_back(urie);
  vriem[el].push_back(vrie);
  vriem[el].push_back(wrie);

  if(e_right != -1)
  {
    vriem[er].push_back(urie);
    vriem[er].push_back(vrie);
    vriem[er].push_back(wrie);
  }
}

std::vector< std::array< tk::real, 3 > >
fluxTerms(
  std::size_t ncomp,
  std::size_t nmat,
  const std::vector< inciter::EOS >& /*mat_blk*/,
  const std::vector< tk::real >& ugp )
// *****************************************************************************
//  Compute the flux-function for the multimaterial PDEs
//! \param[in] ncomp Number of components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
// //! \param[in] mat_blk EOS material block
//! \param[in] ugp State vector at the Gauss point at which flux is required
//! \return Flux vectors for all components in multi-material PDE system
// *****************************************************************************
{
  using inciter::volfracIdx;
  using inciter::densityIdx;
  using inciter::momentumIdx;
  using inciter::energyIdx;
  using inciter::velocityIdx;
  using inciter::pressureIdx;
  using inciter::deformIdx;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  std::vector< std::array< tk::real, 3 > > fl( ncomp, {{0, 0, 0}} );

  auto u = ugp[ncomp+velocityIdx(nmat,0)];
  auto v = ugp[ncomp+velocityIdx(nmat,1)];
  auto w = ugp[ncomp+velocityIdx(nmat,2)];
  if (inciter::haveSolid(nmat, solidx))
  {
    std::vector< tk::real > al(nmat, 0.0);
    std::vector< std::array< std::array< tk::real, 3 >, 3 > > g, asig;
    std::array< std::array< tk::real, 3 >, 3 >
      sig {{ {{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}} }};
    for (std::size_t k=0; k<nmat; ++k)
    {
      al[k] = ugp[volfracIdx(nmat, k)];
      // inv deformation gradient and Cauchy stress tensors
      g.push_back(inciter::getDeformGrad(nmat, k, ugp));
      asig.push_back(inciter::getCauchyStress(nmat, k, ncomp, ugp));
      for (std::size_t i=0; i<3; ++i) asig[k][i][i] -= ugp[ncomp+pressureIdx(nmat,k)];

      for (size_t i=0; i<3; ++i)
        for (size_t j=0; j<3; ++j)
          sig[i][j] += asig[k][i][j];
    }

    // conservative part of momentum flux
    fl[momentumIdx(nmat, 0)][0] = ugp[momentumIdx(nmat, 0)] * u - sig[0][0];
    fl[momentumIdx(nmat, 1)][0] = ugp[momentumIdx(nmat, 1)] * u - sig[0][1];
    fl[momentumIdx(nmat, 2)][0] = ugp[momentumIdx(nmat, 2)] * u - sig[0][2];

    fl[momentumIdx(nmat, 0)][1] = ugp[momentumIdx(nmat, 0)] * v - sig[1][0];
    fl[momentumIdx(nmat, 1)][1] = ugp[momentumIdx(nmat, 1)] * v - sig[1][1];
    fl[momentumIdx(nmat, 2)][1] = ugp[momentumIdx(nmat, 2)] * v - sig[1][2];

    fl[momentumIdx(nmat, 0)][2] = ugp[momentumIdx(nmat, 0)] * w - sig[2][0];
    fl[momentumIdx(nmat, 1)][2] = ugp[momentumIdx(nmat, 1)] * w - sig[2][1];
    fl[momentumIdx(nmat, 2)][2] = ugp[momentumIdx(nmat, 2)] * w - sig[2][2];

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      fl[volfracIdx(nmat, k)][0] = 0.0;
      fl[volfracIdx(nmat, k)][1] = 0.0;
      fl[volfracIdx(nmat, k)][2] = 0.0;

      // conservative part of material continuity flux
      fl[densityIdx(nmat, k)][0] = u * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][1] = v * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][2] = w * ugp[densityIdx(nmat, k)];

      // conservative part of material total-energy flux
      fl[energyIdx(nmat, k)][0] = u * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][0] - v * asig[k][1][0] - w * asig[k][2][0];
      fl[energyIdx(nmat, k)][1] = v * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][1] - v * asig[k][1][1] - w * asig[k][2][1];
      fl[energyIdx(nmat, k)][2] = w * ugp[energyIdx(nmat, k)]
        - u * asig[k][0][2] - v * asig[k][1][2] - w * asig[k][2][2];

      // conservative part of material inverse deformation gradient
      // g_ij: \partial (g_il u_l) / \partial (x_j)
      if (solidx[k] > 0)
      {
        for (std::size_t i=0; i<3; ++i)
        {
          for (std::size_t j=0; j<3; ++j)
          {
            fl[deformIdx(nmat,solidx[k],i,j)][j] =
              u*g[k][i][0] + v*g[k][i][1] + w*g[k][i][2];
          }
          // other components are zero
        }
      }
    }
  }
  else
  {
    tk::real p(0.0);
    std::vector< tk::real > apk( nmat, 0.0 );
    for (std::size_t k=0; k<nmat; ++k)
    {
      apk[k] = ugp[ncomp+pressureIdx(nmat,k)];
      p += apk[k];
    }

    // conservative part of momentum flux
    fl[momentumIdx(nmat, 0)][0] = ugp[momentumIdx(nmat, 0)] * u + p;
    fl[momentumIdx(nmat, 1)][0] = ugp[momentumIdx(nmat, 1)] * u;
    fl[momentumIdx(nmat, 2)][0] = ugp[momentumIdx(nmat, 2)] * u;

    fl[momentumIdx(nmat, 0)][1] = ugp[momentumIdx(nmat, 0)] * v;
    fl[momentumIdx(nmat, 1)][1] = ugp[momentumIdx(nmat, 1)] * v + p;
    fl[momentumIdx(nmat, 2)][1] = ugp[momentumIdx(nmat, 2)] * v;

    fl[momentumIdx(nmat, 0)][2] = ugp[momentumIdx(nmat, 0)] * w;
    fl[momentumIdx(nmat, 1)][2] = ugp[momentumIdx(nmat, 1)] * w;
    fl[momentumIdx(nmat, 2)][2] = ugp[momentumIdx(nmat, 2)] * w + p;

    for (std::size_t k=0; k<nmat; ++k)
    {
      // conservative part of volume-fraction flux
      fl[volfracIdx(nmat, k)][0] = 0.0;
      fl[volfracIdx(nmat, k)][1] = 0.0;
      fl[volfracIdx(nmat, k)][2] = 0.0;

      // conservative part of material continuity flux
      fl[densityIdx(nmat, k)][0] = u * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][1] = v * ugp[densityIdx(nmat, k)];
      fl[densityIdx(nmat, k)][2] = w * ugp[densityIdx(nmat, k)];

      // conservative part of material total-energy flux
      auto hmat = ugp[energyIdx(nmat, k)] + apk[k];
      fl[energyIdx(nmat, k)][0] = u * hmat;
      fl[energyIdx(nmat, k)][1] = v * hmat;
      fl[energyIdx(nmat, k)][2] = w * hmat;
    }
  }

  return fl;
}
}// tk::
