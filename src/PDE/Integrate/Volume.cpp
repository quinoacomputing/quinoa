// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Volume.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing volume integrals for a system of PDEs in DG
     methods
  \details   This file contains functionality for computing volume integrals for
     a system of PDEs used in discontinuous Galerkin methods for various orders
     of numerical representation.
*/
// *****************************************************************************

#include "Volume.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "ViscousTerms.hpp"
#include "MultiMatTerms.hpp"
#include "Kokkos_Core.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

using range_policy = Kokkos::RangePolicy<execution_space>;
using UnManagedMem = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

template <typename T>
auto changeToView(T* object, size_t n) {
    Kokkos::View<T*, Kokkos::LayoutLeft, Kokkos::HostSpace, UnManagedMem> object_view(object, n);
    return object_view;
}

namespace inciter {
  extern ctr::InputDeck g_inputdeck;
};

// Forward declaration of device function
namespace tk {
  KOKKOS_INLINE_FUNCTION
  void update_rhs_device( ncomp_t ncomp,
                  const std::size_t ndof,
                  const std::size_t ndof_el,
                  const tk::real wt,
                  const std::size_t m_nprop,
                  const std::size_t e,
                  const Kokkos::Array<Kokkos::Array<tk::real, NDOF_MAX>, 3>& dBdx,
                  const Kokkos::Array<Kokkos::Array<tk::real, 3>, NCOMP_MAX>& fl,
                  Kokkos::View<real*, memory_space> R);
}

template< class ViscousTerms >
void
tk::volIntViscous(
            const ViscousTerms& viscousRhs,
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
            Fields& R)
// *****************************************************************************
//  Compute volume integrals for viscous DG
//! \tparam ViscousTerms Policy type that computes PDE-specific viscous RHS
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > state(ncomp+nprim);

  Assert( ndof*ncomp == R.nprop(),
          "Mismatch in viscous RHS polynomial and component sizes" );

  std::vector<std::array<tk::real, 3>> visc_fl(ncomp);
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

    auto dof_el = ndofel[e];

    // Compute the derivatives of basis function for second order terms
    std::array< std::vector<tk::real>, 3 > dBdx;
    for (std::size_t i=0; i<3; ++i) dBdx[i].resize( dof_el, 0 );
    eval_dBdx_p1( dof_el, jacInv, dBdx );

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (dof_el > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp],
        B );

      auto wt = wgp[igp] * geoElem(e, 0);

      // volume fluxes
      if(dof_el > 1)
      {
        state = viscousRhs.stateAt(mat_blk, U, P, e, dof_el, B );
        // compute viscous flux
        viscousRhs.volumeFlux( mat_blk, ncomp, state, visc_fl ); //currently no-op

        update_rhs( ncomp, ndof, dof_el, wt, e, dBdx, visc_fl, R );
      }
    }
  }
}

void
tk::volInt( std::size_t nmat,
            real t,
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::size_t rdof,
            const std::size_t nelem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const Fields& geoElem,
            const FluxFn& flux,
            const VelFn& vel,
            const SrcFn& src,
            const Fields& U,
            const Fields& P,
            const Fields& W,
            const std::vector< std::size_t >& ndofel,
            Fields& R,
            int intsharp )
// *****************************************************************************
//  Compute volume integrals for DG
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] src Source function to use
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& ale = inciter::g_inputdeck.get< tag::ale, tag::ale >();
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

    auto dof_el = ndofel[e];

    // Compute the derivatives of basis function for second order terms
    std::array< std::vector<tk::real>, 3 > dBdx;
    for (std::size_t i=0; i<3; ++i) dBdx[i].resize( dof_el, 0 );
    eval_dBdx_p1( dof_el, jacInv, dBdx );

    std::vector< tk::real > B(dof_el);

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      if (dof_el > 4)
        eval_dBdx_p2( igp, coordgp, jacInv, dBdx );

      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      eval_basis( dof_el, coordgp[0][igp], coordgp[1][igp], coordgp[2][igp],
        B );

      auto wt = wgp[igp] * geoElem(e, 0);

      // volume fluxes
      if(dof_el > 1)
      {
        evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, e, ndofel[e], inpoel, coord, geoElem,
          {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U, P, state);

        // evaluate prescribed velocity (if any)
        auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

        // compute flux
        auto fl = flux( ncomp, mat_blk, state, v );

        // update flux according to mesh velocity at quadrature point
        if (ale) {
          auto w_igp = evaluateMeshVelTet( e, igp, inpoel, coordgp, W );

          for (std::size_t c=0; c<ncomp; ++c) {
            for (std::size_t i=0; i<3; ++i) {
              fl[c][i] -= state[c]*w_igp[i];
            }
          }
        }

        update_rhs( ncomp, ndof, dof_el, wt, e, dBdx, fl, R );
      }

      // source terms
      std::vector< real > s(ncomp, 0.0);
      src( nmat, mat_blk, gp[0], gp[1], gp[2], t, s );

      update_rhs_src( ndof, ndofel[e], wt, e, B, s, R );
    }
  }
}

void tk::update_rhs( ncomp_t ncomp,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t e,
                const std::array< std::vector<tk::real>, 3 >& dBdx,
                const std::vector< std::array< tk::real, 3 > >& fl,
                Fields& R )
// *****************************************************************************
//  Update the rhs by adding the flux term integrals (CPU version)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] dBdx Vector of basis function derivatives
//! \param[in] fl Vector of numerical flux
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( dBdx[0].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[1].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( dBdx[2].size() == ndof_el,
    "Size mismatch for basis function derivatives" );
  Assert( fl.size() == ncomp, "Size mismatch for flux term" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e, mark+1) +=
      wt * (fl[c][0]*dBdx[0][1] + fl[c][1]*dBdx[1][1] + fl[c][2]*dBdx[2][1]);
    R(e, mark+2) +=
      wt * (fl[c][0]*dBdx[0][2] + fl[c][1]*dBdx[1][2] + fl[c][2]*dBdx[2][2]);
    R(e, mark+3) +=
      wt * (fl[c][0]*dBdx[0][3] + fl[c][1]*dBdx[1][3] + fl[c][2]*dBdx[2][3]);

    if( ndof_el > 4 )
    {
      R(e, mark+4) +=
        wt * (fl[c][0]*dBdx[0][4] + fl[c][1]*dBdx[1][4] + fl[c][2]*dBdx[2][4]);
      R(e, mark+5) +=
        wt * (fl[c][0]*dBdx[0][5] + fl[c][1]*dBdx[1][5] + fl[c][2]*dBdx[2][5]);
      R(e, mark+6) +=
        wt * (fl[c][0]*dBdx[0][6] + fl[c][1]*dBdx[1][6] + fl[c][2]*dBdx[2][6]);
      R(e, mark+7) +=
        wt * (fl[c][0]*dBdx[0][7] + fl[c][1]*dBdx[1][7] + fl[c][2]*dBdx[2][7]);
      R(e, mark+8) +=
        wt * (fl[c][0]*dBdx[0][8] + fl[c][1]*dBdx[1][8] + fl[c][2]*dBdx[2][8]);
      R(e, mark+9) +=
        wt * (fl[c][0]*dBdx[0][9] + fl[c][1]*dBdx[1][9] + fl[c][2]*dBdx[2][9]);
    }
  }
}

//! Kokkos device version of update_rhs
KOKKOS_INLINE_FUNCTION
void tk::update_rhs_device( ncomp_t ncomp,
                const std::size_t ndof,
                const std::size_t ndof_el,
                const tk::real wt,
                const std::size_t m_nprop,
                const std::size_t e,
                const Kokkos::Array<Kokkos::Array<tk::real, NDOF_MAX>, 3>& dBdx,
                const Kokkos::Array<Kokkos::Array<tk::real, 3>, NCOMP_MAX>& fl,
                Kokkos::View<real*, memory_space> R)
// *****************************************************************************
//  Update the rhs by adding the flux term integrals (Kokkos device version)
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] m_nprop Number of properties in R
//! \param[in] e Element index
//! \param[in] dBdx Array of basis function derivatives
//! \param[in] fl Array of numerical flux
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(e * m_nprop + mark+1) +=
      wt * (fl[c][0]*dBdx[0][1] + fl[c][1]*dBdx[1][1] + fl[c][2]*dBdx[2][1]);
    R(e * m_nprop + mark+2) +=
      wt * (fl[c][0]*dBdx[0][2] + fl[c][1]*dBdx[1][2] + fl[c][2]*dBdx[2][2]);
    R(e * m_nprop + mark+3) +=
      wt * (fl[c][0]*dBdx[0][3] + fl[c][1]*dBdx[1][3] + fl[c][2]*dBdx[2][3]);

    if( ndof_el > 4 )
    {
      R(e * m_nprop + mark+4) +=
        wt * (fl[c][0]*dBdx[0][4] + fl[c][1]*dBdx[1][4] + fl[c][2]*dBdx[2][4]);
      R(e * m_nprop + mark+5) +=
        wt * (fl[c][0]*dBdx[0][5] + fl[c][1]*dBdx[1][5] + fl[c][2]*dBdx[2][5]);
      R(e * m_nprop + mark+6) +=
        wt * (fl[c][0]*dBdx[0][6] + fl[c][1]*dBdx[1][6] + fl[c][2]*dBdx[2][6]);
      R(e * m_nprop + mark+7) +=
        wt * (fl[c][0]*dBdx[0][7] + fl[c][1]*dBdx[1][7] + fl[c][2]*dBdx[2][7]);
      R(e * m_nprop + mark+8) +=
        wt * (fl[c][0]*dBdx[0][8] + fl[c][1]*dBdx[1][8] + fl[c][2]*dBdx[2][8]);
      R(e * m_nprop + mark+9) +=
        wt * (fl[c][0]*dBdx[0][9] + fl[c][1]*dBdx[1][9] + fl[c][2]*dBdx[2][9]);
    }
  }
}

void
tk::update_rhs_src(
  const std::size_t ndof,
  const std::size_t ndof_el,
  const tk::real wt,
  const std::size_t e,
  const std::vector< tk::real >& B,
  const std::vector< tk::real >& s,
  Fields& R )
// *****************************************************************************
//  Update the rhs by adding the source term integrals
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_el Number of degrees of freedom for local element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] e Element index
//! \param[in] B Vector of basis functions
//! \param[in] s Source term vector
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  for (ncomp_t c=0; c<s.size(); ++c)
  {
    auto mark = c*ndof;
    R(e, mark)   += wt * s[c];

    if ( ndof_el > 1 )
    {
      R(e, mark+1) += wt * s[c] * B[1];
      R(e, mark+2) += wt * s[c] * B[2];
      R(e, mark+3) += wt * s[c] * B[3];

      if( ndof_el > 4 )
      {
        R(e, mark+4) += wt * s[c] * B[4];
        R(e, mark+5) += wt * s[c] * B[5];
        R(e, mark+6) += wt * s[c] * B[6];
        R(e, mark+7) += wt * s[c] * B[7];
        R(e, mark+8) += wt * s[c] * B[8];
        R(e, mark+9) += wt * s[c] * B[9];
      }    }
  }
}

void
tk::volInt_constP(
  std::size_t nmat,
  real t,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::size_t nelem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const Fields& geoElem,
  const FluxFn& flux,
  const VelFn& vel,
  const SrcFn& src,
  const Fields& U,
  const Fields& P,
  const Fields& W,
  Fields& R,
  int intsharp,
  VolIntDeviceViews* dev, bool prestaged //added
)
// *****************************************************************************
//  Compute volume integrals for const-order DG with Kokkos acceleration
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Total number of degrees of freedom included reconstructed ones
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] flux Flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] src Source function to use
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
//! \param[in,out] R Right-hand side vector added to
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  // TODO: Wire up the W variable to the current datastructure
  Kokkos::Profiling::pushRegion("volInt");

  const auto& ale = inciter::g_inputdeck.get< tag::ale, tag::ale >();
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  size_t ncomp = U.nprop()/rdof;
  size_t nprim = P.nprop()/rdof;

  size_t m_nprop = U.nprop();
  size_t p_nprop = P.nprop();
  size_t r_nprop = R.nprop();
  size_t geo_nprop = geoElem.nprop();

  const auto& solidx = inciter::g_inputdeck.get< //added the &
      tag::matidxmap, tag::solidx >();

  auto bparam = inciter::g_inputdeck.get< tag::multimat,
    tag::intsharp_param >();

  // Fail loudly on host side if this config overruns the fixed-size device scratch buffers
  // Avoid corruption of private thread memory
  checkKokkosCaps( nmat,ndof,rdof,ncomp,nprim );

  // Quadrature points computed once (constant P)
  auto ng = tk::NGvol(ndof);

  // Persistent device buffers
  // Does not get rid of host->device deep_copy since the data changes per step
  // But removes need for constant cudaMalloc and cudaFree
  VolIntDeviceViews local_dev;
  VolIntDeviceViews& dv = dev ? *dev : local_dev;

  // Check if mesh data currently on the device belongs to this partition at this mesh gen state
  // PDE object lives in global g_dgpde vector and shared by every DG char on the PE
  // Returns false if time-invariant views need to be re-uploaded
  const bool mesh_ok = meshResident( dv, inpoel, coord, geoElem, nelem, nmat );
  
  // Transfer solidx vector
  auto solidx_h_view = changeToView(solidx.data(), nmat);
  if (ensureDeviceCapacity(dv.solidx, "solidx_d_view", nmat) || !mesh_ok)
    Kokkos::deep_copy(dv.solidx, solidx_h_view);

  // Transfer inpoel variable
  size_t inpoel_size = inpoel.size();
  auto inpoel_h_view = changeToView(inpoel.data(), inpoel_size);
  if (ensureDeviceCapacity(dv.inpoel, "inpoel_d_view", inpoel_size) || !mesh_ok)
    Kokkos::deep_copy(dv.inpoel, inpoel_h_view);

  // Transfer coord (nodal coordinates)
  size_t coordx_size = coord[0].size();
  auto cx_h_view = changeToView(coord[0].data(), coordx_size);
  if (ensureDeviceCapacity(dv.cx, "cx_d_view", coordx_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cx, cx_h_view);

  size_t coordy_size = coord[1].size();
  auto cy_h_view = changeToView(coord[1].data(), coordy_size);
  if (ensureDeviceCapacity(dv.cy, "cy_d_view", coordy_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cy, cy_h_view);
  
  size_t coordz_size = coord[2].size();
  auto cz_h_view = changeToView(coord[2].data(), coordz_size);
  if (ensureDeviceCapacity(dv.cz, "cz_d_view", coordz_size) || !mesh_ok)
    Kokkos::deep_copy(dv.cz, cz_h_view);

  // geoElem, U, P, R transfer
  size_t geoElem_size = geoElem.getSize();
  auto geoElem_h_view = changeToView(geoElem.getPointer(), geoElem_size);
  if (ensureDeviceCapacity(dv.geoElem, "geoElem_d_view", geoElem_size) || !mesh_ok)
    Kokkos::deep_copy(dv.geoElem, geoElem_h_view);

  // U,P,R change every call, but when prestaged caller has queued them onto dv
  // on the same exec space instance, while also owning the R round trip
  // hence, we can save on an upload/download as R arrives carrying surf int contributions
  auto exec = Kokkos::DefaultExecutionSpace();
  const size_t R_size = R.getSize();

  if(!prestaged){
    const size_t P_size = P.getSize();
    uploadStaged( exec, dv.P, dv.stage_P, P.getPointer(), P_size, "P_d_view" );
    
    const size_t U_size = U.getSize();
    uploadStaged( exec, dv.U, dv.stage_U, U.getPointer(), U_size, "U_d_view" );
#ifdef VOLINT_R_PRESEEDED
    uploadStaged( exec, dv.R, dv.stage_R, R.getPointer(), R_size, "R_d_view" );
#else
    ensureDeviceCapacity(dv.R, "R_d_view", R_size);
    Kokkos::deep_copy(exec, dv.R, 0.0);
#endif
  }

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
  auto R_d_view = dv.R;

  // Quadrature points
  // Can be hoisted out because P is constant
  Kokkos::Array<Kokkos::Array<real, NQUAD_MAX>, 3> coordgp = {};
  Kokkos::Array<real, NQUAD_MAX> wgp = {};
  GaussQuadratureTet(ng, coordgp, wgp );

  Kokkos::parallel_for("volInt_kernel",range_policy(0, nelem), KOKKOS_LAMBDA(const size_t e)
  {
    if(ndof > 1)
    {
      // Extract the element coordinates
      Kokkos::Array<Kokkos::Array<real, 3>, 4> coordel;
      for (int i=0; i<4; i++) {
        coordel[i][0] = cx_d_view(inpoel_d_view(4*e + i));
        coordel[i][1] = cy_d_view(inpoel_d_view(4*e + i));
        coordel[i][2] = cz_d_view(inpoel_d_view(4*e + i));
      }

      Kokkos::Array<real, NDOF_MAX> B = {};
      Kokkos::Array<Kokkos::Array<real, NDOF_MAX>, 3> dBdx = {};
      auto jacInv =
        inverseJacobian(coordel[0], coordel[1], coordel[2], coordel[3] );
 
      eval_dBdx_p1(ndof, jacInv, dBdx);

      Kokkos::Array<Kokkos::Array<Kokkos::Array<real, 3>, 3>, NMAT_MAX> g = {};
      Kokkos::Array<Kokkos::Array<Kokkos::Array<real, 3>, 3>, NMAT_MAX> asig = {};
      Kokkos::Array<real, NMAT_MAX> al = {};
      Kokkos::Array<Kokkos::Array<real, 3>, NCOMP_MAX> fl = {};
      Kokkos::Array<real, NMAT_MAX> apk = {};
      Kokkos::Array<real, NSTATE_MAX> state = {};

      for (std::size_t igp=0; igp<ng; ++igp)
      {
        if (ndof > 4)
          eval_dBdx_p2( igp, coordgp, jacInv, dBdx);

        // Compute the coordinates of quadrature point at physical domain
        // auto gp = eval_gp( igp, coordel, coordgp);
  
        // Compute the basis function
        eval_basis( rdof, coordgp[0][igp], coordgp[1][igp],
                    coordgp[2][igp], B);

        auto wt = wgp[igp] * geoElem_d_view(e * geo_nprop);
        
        evalPolynomialSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, e, rdof, m_nprop, p_nprop, geo_nprop,
          bparam, solidx_d_view, inpoel_d_view, 
          cx_d_view, cy_d_view, cz_d_view, geoElem_d_view,
          {{coordgp[0][igp], coordgp[1][igp], coordgp[2][igp]}}, B, U_d_view, 
          P_d_view, state);

        // compute flux
        fluxTerms_multimat_kokkos(ncomp, nmat, solidx_d_view, 
          mat_blk, state, g, asig, al, fl, apk);
            
        // Call device version explicitly
        tk::update_rhs_device(ncomp, ndof, ndof, wt, r_nprop, e, dBdx, fl, R_d_view);  
      }
    }
  });
  
  //Kokkos::deep_copy(R_h_view, R_d_view);

  // Queued on same instance as kernel so ordered after it
  // downloadStaged() fences before copying out of pinned buffer
  // and we skip when prestaged because the caller downloads R once after last kernel
  if(!prestaged)
    downloadStaged( exec, R.getPointerNonConst(), dv.stage_R, R_d_view, R_size, "R_d_view" );

  Kokkos::Profiling::popRegion();
}

void
tk::srcInt_constP( std::size_t nmat,
                   real t,
                   const std::vector< inciter::EOS >& mat_blk,
                   const std::size_t ndof,
                   const std::size_t ncomp,
                   const std::size_t nelem,
                   const std::vector< std::size_t >& inpoel,
                   const UnsMesh::Coords& coord,
                   const Fields& geoElem,
                   const SrcFn& src,
                   Fields& R)
//  ****************************************************************************
//  Compute source terms integrals for const-order DG (not p-adaptive)
//! \details Split out of volInt_constP(). src() is a host std::function and is
//!   not device-callable so this contribution is applied on the host.
//!   Pulled out separately to make explicit that rhs() must call this AFTER 
//!   the last D2H of R, which will make sense once rhs() manages R across kernels.
//  ****************************************************************************
{
  Kokkos::Profiling::pushRegion("srcInt_constP");

  const auto ng = tk::NGvol(ndof);

  std::array< std::vector<real>,3 > coordgp_h;
  std::vector<real> wgp_h;
  for (std::size_t i=0; i<3; ++i) coordgp_h[i].resize(ng);
  wgp_h.resize(ng);
  GaussQuadratureTet(ng,coordgp_h,wgp_h);

  // Dubiner basis is a function of ref coords only, so constP implies identical for all elems
  std::vector<std::vector<real>> Bg(ng,std::vector<real>(ndof));
  for (std::size_t igp=0; igp<ng; ++igp)
    eval_basis(ndof, coordgp_h[0][igp], coordgp_h[1][igp], coordgp_h[2][igp], Bg[igp]);

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  std::vector<real> sv(ncomp,0.0);
    
  for (std::size_t e=0; e<nelem; ++e)
  {
    std::array< std::array<real,3>, 4 > coordel {{
      {{ cx[inpoel[4*e]], cy[inpoel[4*e]], cz[inpoel[4*e]] }},
      {{ cx[inpoel[4*e+1]], cy[inpoel[4*e+1]], cz[inpoel[4*e+1]] }},
      {{ cx[inpoel[4*e+2]], cy[inpoel[4*e+2]], cz[inpoel[4*e+2]] }},
      {{ cx[inpoel[4*e+3]], cy[inpoel[4*e+3]], cz[inpoel[4*e+3]] }}
    }};

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      auto gp = eval_gp( igp, coordel, coordgp_h );
      auto wt = wgp_h[igp]*geoElem(e,0);
      
      std::fill( begin(sv), end(sv), 0.0 );
      src( nmat, mat_blk, gp[0], gp[1], gp[2], t, sv );
      update_rhs_src( ndof, ndof, wt, e, Bg[igp], sv, R );
    }
  }

  Kokkos::Profiling::popRegion();
}

void tk::srcIntFV( const std::vector< inciter::EOS >& mat_blk,
              real t,
              const std::size_t nelem,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R,
              std::size_t nmat )
// *****************************************************************************
//  Compute source term integrals for DG
//! \param[in] mat_blk Material block EOS
//! \param[in] t Physical time
//! \param[in] nelem Maximum number of elements
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in,out] R Right-hand side vector computed
//! \param[in] nmat Number of materials. A default is set to 1, so that calling
//!   code for single material systems primitive quantities does not need to
//!   specify this argument.
// *****************************************************************************
{
  auto ncomp = R.nprop();

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Compute the source term variable
    std::vector< real > s(ncomp, 0.0);
    src( nmat, mat_blk, geoElem(e,1), geoElem(e,2), geoElem(e,3), t, s );

    // Add the source term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(e, c) += geoElem(e,0) * s[c];
    }
  }
}

void
tk::volIntViscousMultiSpecies(
  std::size_t nspec,
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
  Fields& R )
// *****************************************************************************
//  Compute volume integrals of viscous fluxes for multispecies flow
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  if (ndof == 4) {
    MultiSpeciesViscousTermsDGP1 viscousRhs( nspec, rdof );
    volIntViscous( viscousRhs, mat_blk, ndof, rdof, nelem,
      inpoel, coord, geoElem, U, P, ndofel, R );
  }
  else if (ndof ==1) {
     //Do nothing but don't exit as nothing was here in DGMultiSpecies.cpp before adding this call for DGP1
  }
  else
    Throw( "Viscous operators only implemented for scheme = 'dgp1'." );
}
