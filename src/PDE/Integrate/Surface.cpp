// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing internal surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing internal surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Surface.hpp"
#include "ViscousTerms.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "EoS/GetMatProp.hpp"
#include "Riemann/HLLCMultiMatConstP.hpp"
#include "Inciter/Options/Flux.hpp"

namespace {
// Accumulators used by update_rhs_fa_constP so that one templated body can write either
// to host tk::Fields/std::vector or to device Kokkos::Views
// Device overloads are atomic so an internal face is written by both of its adjacent elems
  
using UnManagedMem = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

// Wrap a host raw pointer in an unmanaged host view
// Mirrors the helper of the same name in MultiMatTerms.cpp, which lives in an
// anonymous namespace there and so is not reachable from here.
template< typename T >
auto changeToView( T* object, std::size_t n ) {
  Kokkos::View< T*, Kokkos::LayoutLeft, Kokkos::HostSpace, UnManagedMem >
    object_view( object, n );
  return object_view;
}

// Accumulate into R, host tk::Fields
inline void
rhsAccum( tk::Fields& R, std::size_t e, std::size_t /*nprop*/,
          std::size_t idx, tk::real v )
{
  R(e,idx) += v;
}  

// Accumulate into R, Kokkos::Views (device)
KOKKOS_INLINE_FUNCTION
void
rhsAccum( const Kokkos::View< tk::real*, memory_space >& R, std::size_t e,
          std::size_t nprop, std::size_t idx, tk::real v )
{
  Kokkos::atomic_add( &R(e*nprop + idx), v );
}

// Accumulate into riemannDeriv, host vector of vectors
inline void
rdAccum( std::vector< std::vector< tk::real > >& rd, std::size_t row,
         std::size_t col, std::size_t /*ncol*/, tk::real v )
{
  rd[row][col] += v;
}

// Accumulate into riemannDeriv, device view (atomic)
KOKKOS_INLINE_FUNCTION
void
rdAccum( const Kokkos::View< tk::real*, memory_space >& rd, std::size_t row,
         std::size_t col, std::size_t ncol, tk::real v )
{
  Kokkos::atomic_add( &rd(row*ncol + col), v );
}
}//namespace

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace tk {

namespace {

template< class ViscousTerms >
void
viscousInternalFaceInt(
  const ViscousTerms& viscousRhs,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  Fields& R )
// *****************************************************************************
//! \brief Shared PDE-nonspecific face traversal for viscous surface operators
//! \tparam ViscousTerms Policy type that computes PDE-specific viscous RHS
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Number of active solution degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in,out] R Right-hand side vector computed
//! \details This routine owns only the face traversal and local geometric
//!   quantities common to viscous surface terms. The supplied policy receives
//!   face-local data and provides the PDE-specific flux.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];
  const auto ncomp = static_cast< ncomp_t >( R.nprop()/ndof );

  Assert( ndof*ncomp == R.nprop(),
          "Mismatch in viscous RHS polynomial and component sizes" );

  std::array< std::vector< tk::real >, 2 > B;
  std::array< std::vector< tk::real >, 2 > Bcc;
  std::array< std::vector< tk::real >, 2 > state;
  std::array< std::vector< tk::real >, 2 > cellAvgState;
  std::array< std::array< std::vector< tk::real >, 3 >, 2 > dBdx;
  std::array< std::array< std::array< tk::real, 3 >, 4>, 2 > grad;
  std::vector< tk::real > fl( ncomp, 0.0 );

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);
    auto ndof_l = viscousRhs.localDof( el );
    auto ndof_r = viscousRhs.localDof( er );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // face normal
    std::array< real, 3 > fn{{geoFace(f,1), geoFace(f,2), geoFace(f,3)}};

    // face centroid
    std::array< real, 3 > gp{{geoFace(f,4), geoFace(f,5), geoFace(f,6)}};

    // Quadrature points in element-reference-frames
    std::array< tk::real, 3> ref_gp_l{{
      Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
      Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l }};
    std::array< tk::real, 3> ref_gp_r{{
      Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
      Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r }};

    // Compute the basis functions
    B[0].resize(ndof_l);
    B[1].resize(ndof_r);
    eval_basis( ndof_l, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B[0] );
    eval_basis( ndof_r, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2], B[1] );

    std::array< std::array< real, 3 >, 2 > centroids{{
      {{geoElem(el,1), geoElem(el,2), geoElem(el,3)}},
      {{geoElem(er,1), geoElem(er,2), geoElem(er,3)}} }};

    // Gradients of basis functions
    for (std::size_t i=0; i<3; ++i) {
      dBdx[0][i].assign( ndof_l, 0.0 );
      dBdx[1][i].assign( ndof_r, 0.0 );
    }
    auto jacInv_l =
      inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    eval_dBdx_p1( ndof_l, jacInv_l, dBdx[0] );
    auto jacInv_r =
      inverseJacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );
    eval_dBdx_p1( ndof_r, jacInv_r, dBdx[1] );

    // Compute high-order face states
    state[0] = viscousRhs.stateAt( mat_blk, U, P, el, ndof_l, B[0] );
    state[1] = viscousRhs.stateAt( mat_blk, U, P, er, ndof_r, B[1] );

    // Compute cell-average states
    Bcc[0].assign( ndof_l, 0.0 );
    Bcc[1].assign( ndof_r, 0.0 );
    Bcc[0][0] = 1.0;
    Bcc[1][0] = 1.0;
    cellAvgState[0] =
      viscousRhs.stateAt( mat_blk, U, P, el, ndof_l, Bcc[0] );
    cellAvgState[1] =
      viscousRhs.stateAt( mat_blk, U, P, er, ndof_r, Bcc[1] );

    // Compute gradients
    viscousRhs.gradientIntElem( U, P, el, dBdx[0], grad[0] );
    viscousRhs.gradientIntElem( U, P, er, dBdx[1], grad[1] );

    // Compute viscous fluxes
    viscousRhs.interiorFlux( mat_blk, ncomp, state, cellAvgState, fn,
      centroids, grad, fl );

    // Contribute fluxes to RHS
    for (ncomp_t c=0; c<ncomp; ++c) {
      R(el, c) += geoFace(f,0) * fl[c];
      R(er, c) -= geoFace(f,0) * fl[c];
    }
  }
}

template< class ViscousTerms >
void
viscousInternalFaceIntDG(
  const ViscousTerms& viscousRhs,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  Fields& R  )
// *****************************************************************************
//  Compute internal surface flux integrals
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{

  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  const auto ncomp = static_cast< ncomp_t >( R.nprop()/ndof );

  Assert( ndof*ncomp == R.nprop(),
          "Mismatch in viscous RHS polynomial and component sizes" );

  std::array< std::vector< tk::real >, 2 > B;
  std::array< std::vector< tk::real >, 2 > state;
  std::array< std::array< std::vector< tk::real >, 3 >, 2 > dBdx;
  std::array< std::array< std::array< tk::real, 3 >, 4>, 2 > grad;
  std::vector< tk::real > fl( ncomp, 0.0 );

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);
    auto ndof_l = viscousRhs.localDof( el );
    auto ndof_r = viscousRhs.localDof( er );

    auto ng_l = tk::NGfa(ndof_l);
    auto ng_r = tk::NGfa(ndof_r);

    // When the number of gauss points for the left and right element are
    // different, choose the larger ng
    auto ng = std::max( ng_l, ng_r );

    // arrays for quadrature points
    std::array< std::vector< real >, 2 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    // get quadrature point weights and coordinates for triangle
    GaussQuadratureTri( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // Extract the face coordinates
    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
      {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
      {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }} }};

    std::array< real, 3 >
      fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordfa, coordgp );

      std::array< tk::real, 3> ref_gp_l{
        Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };
      std::array< tk::real, 3> ref_gp_r{
        Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r };

      // Compute the basis functions
      B[0].resize(ndof_l);
      B[1].resize(ndof_r);
      auto B_l = B[0];
      auto B_r = B[1];

      eval_basis( ndof_l, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );
      eval_basis( ndof_r, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2], B_r );

      auto wt = wgp[igp] * geoFace(f,0);

      // Gradients of basis functions
      for (std::size_t i=0; i<3; ++i) {
        dBdx[0][i].assign( ndof_l, 0.0 );
        dBdx[1][i].assign( ndof_r, 0.0 );
      }

      auto jacInv_l =
        inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
      eval_dBdx_p1( ndof_l, jacInv_l, dBdx[0] );
      auto jacInv_r =
        inverseJacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );
      eval_dBdx_p1( ndof_r, jacInv_r, dBdx[1] );

      // Compute high-order face states (reconstruction)
      state[0] = viscousRhs.stateAt( mat_blk, U, P, el, ndof_l, B_l );
      state[1] = viscousRhs.stateAt( mat_blk, U, P, er, ndof_r, B_r );

      // Compute gradients
      viscousRhs.gradientIntElem( U, P, el, dBdx[0], grad[0] ); // currently no-op
      viscousRhs.gradientIntElem( U, P, er, dBdx[1], grad[1] ); // currently no-op

      // Compute viscous fluxes
      viscousRhs.interiorFlux( mat_blk, ncomp, state, fn,
        grad, fl ); // currently no-op

      // Contribute fluxes to RHS
      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        R(el, mark) += wt * fl[c];
        R(er, mark) -= wt * fl[c];

        if(ndof_l > 1) //DG(P1)
        {
          R(el, mark+1) += wt * fl[c] * B_l[1];
          R(el, mark+2) += wt * fl[c] * B_l[2];
          R(el, mark+3) += wt * fl[c] * B_l[3];
        }

        if(ndof_r > 1) //DG(P1)
        {
          R(er, mark+1) -= wt * fl[c] * B_r[1];
          R(er, mark+2) -= wt * fl[c] * B_r[2];
          R(er, mark+3) -= wt * fl[c] * B_r[3];
        }

        if(ndof_l > 4) //DG(P2)
        {
          R(el, mark+4) += wt * fl[c] * B_l[4];
          R(el, mark+5) += wt * fl[c] * B_l[5];
          R(el, mark+6) += wt * fl[c] * B_l[6];
          R(el, mark+7) += wt * fl[c] * B_l[7];
          R(el, mark+8) += wt * fl[c] * B_l[8];
          R(el, mark+9) += wt * fl[c] * B_l[9];
        }

        if(ndof_r > 4) //DG(P2)
        {
          R(er, mark+4) -= wt * fl[c] * B_r[4];
          R(er, mark+5) -= wt * fl[c] * B_r[5];
          R(er, mark+6) -= wt * fl[c] * B_r[6];
          R(er, mark+7) -= wt * fl[c] * B_r[7];
          R(er, mark+8) -= wt * fl[c] * B_r[8];
          R(er, mark+9) -= wt * fl[c] * B_r[9];
        }
      }
    }
  }
}
} // anonymous namespace

void
surfInt( const bool pref,
         std::size_t nmat,
         const std::vector< inciter::EOS >& mat_blk,
         real t,
         const std::size_t ndof,
         const std::size_t rdof,
         const std::vector< std::size_t >& inpoel,
         const std::vector< std::size_t >& /*solidx*/,
         const UnsMesh::Coords& coord,
         const inciter::FaceData& fd,
         const Fields& geoFace,
         const Fields& geoElem,
         const RiemannFluxFn& flux,
         const VelFn& vel,
         const Fields& U,
         const Fields& P,
         const Fields& W,
         const std::vector< std::size_t >& ndofel,
         const tk::real /*dt*/,
         Fields& R,
         std::vector< std::vector< tk::real > >& riemannDeriv,
         int intsharp )
// *****************************************************************************
//  Compute internal surface flux integrals
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] t Physical time
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
// //! \param[in] solidx Material index indicator
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
// //! \param[in] dt Delta time
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& ale = inciter::g_inputdeck.get< tag::ale, tag::ale >();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::array< std::vector< tk::real >, 2 > state;
  state[0].resize(ncomp+nprim);
  state[1].resize(ncomp+nprim);

  //// Determine if we have solids in our problem
  //bool haveSolid = inciter::haveSolid(nmat, solidx);

  //Assert( (nmat==1 ? riemannDeriv.empty() : true), "Non-empty Riemann "
  //        "derivative vector for single material compflow" );

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    auto ng_l = tk::NGfa(ndofel[el]);
    auto ng_r = tk::NGfa(ndofel[er]);

    // When the number of gauss points for the left and right element are
    // different, choose the larger ng
    auto ng = std::max( ng_l, ng_r );

    // arrays for quadrature points
    std::array< std::vector< real >, 2 > coordgp;
    std::vector< real > wgp;

    coordgp[0].resize( ng );
    coordgp[1].resize( ng );
    wgp.resize( ng );

    // get quadrature point weights and coordinates for triangle
    GaussQuadratureTri( ng, coordgp, wgp );

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // Extract the face coordinates
    std::array< std::array< tk::real, 3>, 3 > coordfa {{
      {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
      {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
      {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }} }};

    std::array< real, 3 >
      fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordfa, coordgp );

      // In order to determine the high-order solution from the left and right
      // elements at the surface quadrature points, the basis functions from
      // the left and right elements are needed. For this, a transformation to
      // the reference coordinates is necessary, since the basis functions are
      // defined on the reference tetrahedron only.
      // The transformation relations are shown below:
      //  xi   = Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT
      //  eta  = Jacobian( coordel[0], coordel[2], gp, coordel[3] ) / detT
      //  zeta = Jacobian( coordel[0], coordel[2], coordel[3], gp ) / detT

      // If an rDG method is set up (P0P1), then, currently we compute the P1
      // basis functions and solutions by default. This implies that P0P1 is
      // unsupported in the p-adaptive DG (PDG). This is a workaround until we
      // have rdofel, which is needed to distinguish between ndofs and rdofs per
      // element for pDG.
      std::size_t dof_el, dof_er;
      if (rdof > ndof)
      {
        dof_el = rdof;
        dof_er = rdof;
      }
      else
      {
        dof_el = ndofel[el];
        dof_er = ndofel[er];
      }

      // For multi-material p-adaptive simulation, when dofel = 1, p0p1 is
      // applied and ndof for solution evaluation should be 4
      if(ncomp > 5 && pref) {
        if(dof_el == 1)
          dof_el = 4;
        if(dof_er == 1)
          dof_er = 4;
      }

      std::array< tk::real, 3> ref_gp_l{
        Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };
      std::array< tk::real, 3> ref_gp_r{
        Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r };

      //Compute the basis functions
      std::vector< tk::real > B_l(dof_el), B_r(dof_er);
      eval_basis( dof_el, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );
      eval_basis( dof_er, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2], B_r );

      auto wt = wgp[igp] * geoFace(f,0);

      evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof, nmat, el, dof_el,
        inpoel, coord, geoElem, ref_gp_l, B_l, U, P, state[0]);
      evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof, nmat, er, dof_er,
        inpoel, coord, geoElem, ref_gp_r, B_r, U, P, state[1]);

      // evaluate prescribed velocity (if any)
      auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

      // mesh velocity at quadrature point
      tk::real wn_igp(0.0);
      if (ale) {
        auto w_igp = evaluateMeshVelTri( f, igp, inpofa, coordgp, W );
        // mesh velocity normal to element face
        wn_igp = tk::dot(w_igp, fn);
      }

      // compute flux
      auto fl = flux( mat_blk, fn, state, v, wn_igp );

      // Add the surface integration term to the rhs
      update_rhs_fa( ncomp, nmat, ndof, ndofel[el], ndofel[er], wt, fn,
                     el, er, fl, B_l, B_r, R, riemannDeriv );
    }
  }
}

void
update_rhs_fa( ncomp_t ncomp,
               std::size_t nmat,
               const std::size_t ndof,
               const std::size_t ndof_l,
               const std::size_t ndof_r,
               const tk::real wt,
               const std::array< tk::real, 3 >& fn,
               const std::size_t el,
               const std::size_t er,
               const std::vector< tk::real >& fl,
               const std::vector< tk::real >& B_l,
               const std::vector< tk::real >& B_r,
               Fields& R,
               std::vector< std::vector< tk::real > >& riemannDeriv )
// *****************************************************************************
//  Update the rhs by adding the surface integration term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_l Number of degrees of freedom for left element
//! \param[in] ndof_r Number of degrees of freedom for right element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] fn Face/Surface normal
//! \param[in] el Left element index
//! \param[in] er Right element index
//! \param[in] fl Surface flux
//! \param[in] B_l Basis function for the left element
//! \param[in] B_r Basis function for the right element
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  // following lines commented until rdofel is made available.
  //Assert( B_l.size() == ndof_l, "Size mismatch" );
  //Assert( B_r.size() == ndof_r, "Size mismatch" );

  using inciter::newSolidsAccFn;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(el, mark) -= wt * fl[c];
    R(er, mark) += wt * fl[c];

    if(ndof_l > 1)          //DG(P1)
    {
      R(el, mark+1) -= wt * fl[c] * B_l[1];
      R(el, mark+2) -= wt * fl[c] * B_l[2];
      R(el, mark+3) -= wt * fl[c] * B_l[3];
    }

    if(ndof_r > 1)          //DG(P1)
    {
      R(er, mark+1) += wt * fl[c] * B_r[1];
      R(er, mark+2) += wt * fl[c] * B_r[2];
      R(er, mark+3) += wt * fl[c] * B_r[3];
    }

    if(ndof_l > 4)          //DG(P2)
    {
      R(el, mark+4) -= wt * fl[c] * B_l[4];
      R(el, mark+5) -= wt * fl[c] * B_l[5];
      R(el, mark+6) -= wt * fl[c] * B_l[6];
      R(el, mark+7) -= wt * fl[c] * B_l[7];
      R(el, mark+8) -= wt * fl[c] * B_l[8];
      R(el, mark+9) -= wt * fl[c] * B_l[9];
    }

    if(ndof_r > 4)          //DG(P2)
    {
      R(er, mark+4) += wt * fl[c] * B_r[4];
      R(er, mark+5) += wt * fl[c] * B_r[5];
      R(er, mark+6) += wt * fl[c] * B_r[6];
      R(er, mark+7) += wt * fl[c] * B_r[7];
      R(er, mark+8) += wt * fl[c] * B_r[8];
      R(er, mark+9) += wt * fl[c] * B_r[9];
    }
  }

  // Prep for non-conservative terms in multimat
  if (fl.size() > ncomp)
  {
    // Gradients of partial pressures
    for (std::size_t k=0; k<nmat; ++k)
    {
      for (std::size_t idir=0; idir<3; ++idir)
      {
        riemannDeriv[3*k+idir][el] += wt * fl[ncomp+k] * fn[idir];
        riemannDeriv[3*k+idir][er] -= wt * fl[ncomp+k] * fn[idir];
      }
    }

    // Divergence of velocity multiples basis fucntion( d(uB) / dx )
    for(std::size_t idof = 0; idof < ndof_l; idof++) {
      riemannDeriv[3*nmat+idof][el] += wt * fl[ncomp+nmat] * B_l[idof];
    }
    for(std::size_t idof = 0; idof < ndof_r; idof++) {
      riemannDeriv[3*nmat+idof][er] -= wt * fl[ncomp+nmat] * B_r[idof];
    }

    // Divergence of asigma: d(asig_ij)/dx_j
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
      {
        std::size_t mark = ncomp+nmat+1+3*(solidx[k]-1);

        for (std::size_t i=0; i<3; ++i)
        {
          riemannDeriv[3*nmat+ndof+3*(solidx[k]-1)+i][el] -=
            wt * fl[mark+i];
          riemannDeriv[3*nmat+ndof+3*(solidx[k]-1)+i][er] +=
            wt * fl[mark+i];
        }
      }

    // Derivatives of g: d(g_il)/d(x_j)-d(g_ij)/d(x_l)
    // for i=1,2,3; j=1,2,3; l=1,2,3. Total = 3x3x3 (per solid)
    std::size_t nsld = inciter::numSolids(nmat, solidx);
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t l=0; l<3; ++l)
              if (j != l)
              {
                std::size_t mark1 = ncomp+nmat+1+3*nsld+9*(solidx[k]-1)+3*i+l;
                std::size_t mark2 = ncomp+nmat+1+3*nsld+9*(solidx[k]-1)+3*i+j;
                riemannDeriv[3*nmat+ndof+3*nsld+newSolidsAccFn(k,i,j,l)][el] -=
                  wt * ( fl[mark1] * fn[j] - fl[mark2] * fn[l]);
                riemannDeriv[3*nmat+ndof+3*nsld+newSolidsAccFn(k,i,j,l)][er] +=
                  wt * ( fl[mark1] * fn[j] - fl[mark2] * fn[l]);
              }
  }
}

// *****************************************************************************
//! Update the rhs by adding surface integration term (const-order path)
//! \details Copy of update_rhs_fa, templated on its container types so one body
//!   serves the host loop and the device kernel. Two differences: solidx and
//!   nsld are passed in rather than read from g_inputdeck (a lookup per
//!   quadrature point, and unreachable from device code), and the logical
//!   length of fl arrives as nflx, since a fixed-size buffer's size() is its
//!   capacity. update_rhs_fa itself is untouched, so the p-adaptive path is
//!   unaffected.
// *****************************************************************************
template< class FnT, class FlT, class BT, class SolidxT, class RT, class RiemannDerivT >
KOKKOS_INLINE_FUNCTION
void
update_rhs_fa_constP( ncomp_t ncomp,
               std::size_t nmat,
               const std::size_t ndof,
               const std::size_t ndof_l,
               const std::size_t ndof_r,
               const tk::real wt,
               const FnT& fn,
               const std::size_t el,
               const std::size_t er,
               const FlT& fl,
               const std::size_t nflx,
               const BT& B_l,
               const BT& B_r,
               const SolidxT& solidx,
               const std::size_t nsld,
               RT& R,
               const std::size_t r_nprop,
               RiemannDerivT& riemannDeriv,
               const std::size_t rd_ncol )
// *****************************************************************************
//  Update the rhs by adding the surface integration term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_l Number of degrees of freedom for left element
//! \param[in] ndof_r Number of degrees of freedom for right element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] fn Face/Surface normal
//! \param[in] el Left element index
//! \param[in] er Right element index
//! \param[in] fl Surface flux
//! \param[in] B_l Basis function for the left element
//! \param[in] B_r Basis function for the right element
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  // following lines commented until rdofel is made available.
  //Assert( B_l.size() == ndof_l, "Size mismatch" );
  //Assert( B_r.size() == ndof_r, "Size mismatch" );

  using inciter::newSolidsAccFn;

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    rhsAccum( R, el, r_nprop, mark, -wt*fl[c] );
    rhsAccum( R, er, r_nprop, mark, wt*fl[c] );

    if(ndof_l > 1)          //DG(P1)
    {
      rhsAccum( R, el, r_nprop, mark+1, -wt * fl[c] * B_l[1] );
      rhsAccum( R, el, r_nprop, mark+2, -wt * fl[c] * B_l[2] );
      rhsAccum( R, el, r_nprop, mark+3, -wt * fl[c] * B_l[3] );
    }

    if(ndof_r > 1)          //DG(P1)
    {
      rhsAccum( R, er, r_nprop, mark+1, wt * fl[c] * B_r[1] );
      rhsAccum( R, er, r_nprop, mark+2, wt * fl[c] * B_r[2] );
      rhsAccum( R, er, r_nprop, mark+3, wt * fl[c] * B_r[3] );
    }

    if(ndof_l > 4)          //DG(P2)
    {
      rhsAccum( R, el, r_nprop, mark+4, -wt * fl[c] * B_l[4] );
      rhsAccum( R, el, r_nprop, mark+5, -wt * fl[c] * B_l[5] );
      rhsAccum( R, el, r_nprop, mark+6, -wt * fl[c] * B_l[6] );
      rhsAccum( R, el, r_nprop, mark+7, -wt * fl[c] * B_l[7] );
      rhsAccum( R, el, r_nprop, mark+8, -wt * fl[c] * B_l[8] );
      rhsAccum( R, el, r_nprop, mark+9, -wt * fl[c] * B_l[9] );
    }

    if(ndof_r > 4)          //DG(P2)
    {
      rhsAccum( R, er, r_nprop, mark+4, wt * fl[c] * B_r[4] );
      rhsAccum( R, er, r_nprop, mark+5, wt * fl[c] * B_r[5] );
      rhsAccum( R, er, r_nprop, mark+6, wt * fl[c] * B_r[6] );
      rhsAccum( R, er, r_nprop, mark+7, wt * fl[c] * B_r[7] );
      rhsAccum( R, er, r_nprop, mark+8, wt * fl[c] * B_r[8] );
      rhsAccum( R, er, r_nprop, mark+9, wt * fl[c] * B_r[9] );
    }
  }

  // Prep for non-conservative terms in multimat
  if (nflx > ncomp)
  {
    // Gradients of partial pressures
    for (std::size_t k=0; k<nmat; ++k)
    {
      for (std::size_t idir=0; idir<3; ++idir)
      {
        rdAccum( riemannDeriv, 3*k+idir, el, rd_ncol,
            wt * fl[ncomp+k] * fn[idir] );
        rdAccum( riemannDeriv, 3*k+idir, er, rd_ncol,
            -wt * fl[ncomp+k] * fn[idir] );
      }
    }

    // Divergence of velocity multiples basis fucntion( d(uB) / dx )
    for(std::size_t idof = 0; idof < ndof_l; idof++) {
      rdAccum( riemannDeriv, 3*nmat+idof, el, rd_ncol,
               wt * fl[ncomp+nmat] * B_l[idof] );
    }
    for(std::size_t idof = 0; idof < ndof_r; idof++) {
      rdAccum( riemannDeriv, 3*nmat+idof, er, rd_ncol,
               -wt * fl[ncomp+nmat] * B_r[idof] );
    }

    // Divergence of asigma: d(asig_ij)/dx_j
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
      {
        std::size_t mark = ncomp+nmat+1+3*(solidx[k]-1);

        for (std::size_t i=0; i<3; ++i)
        {
          rdAccum( riemannDeriv, 3*nmat+ndof+3*(solidx[k]-1)+i, el, rd_ncol,
                   -wt * fl[mark+i] );
          rdAccum( riemannDeriv, 3*nmat+ndof+3*(solidx[k]-1)+i, er, rd_ncol,
                   wt * fl[mark+i] );
        }
      }

    // Derivatives of g: d(g_il)/d(x_j)-d(g_ij)/d(x_l)
    // for i=1,2,3; j=1,2,3; l=1,2,3. Total = 3x3x3 (per solid)
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            for (std::size_t l=0; l<3; ++l)
              if (j != l)
              {
                std::size_t mark1 = ncomp+nmat+1+3*nsld+9*(solidx[k]-1)+3*i+l;
                std::size_t mark2 = ncomp+nmat+1+3*nsld+9*(solidx[k]-1)+3*i+j;
                rdAccum( riemannDeriv, 3*nmat+ndof+3*nsld+newSolidsAccFn(k,i,j,l), el, rd_ncol,
                         -wt * ( fl[mark1] * fn[j] - fl[mark2] * fn[l]) );
                rdAccum( riemannDeriv, 3*nmat+ndof+3*nsld+newSolidsAccFn(k,i,j,l), er, rd_ncol,
                         wt * ( fl[mark1] * fn[j] - fl[mark2] * fn[l]) );
              }
  }
}

void
surfInt_constP(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  real t,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& inpoel,
  const std::vector< std::size_t >& solidx,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const RiemannFluxFn& /*flux*/,
  const VelFn& vel,
  const Fields& U,
  const Fields& P,
  const Fields& W,
  const tk::real /*dt*/,
  Fields& R,
  std::vector< std::vector< tk::real > >& riemannDeriv,
  int intsharp,
  [[maybe_unused]] SurfIntDeviceViews* dev,
  [[maybe_unused]] bool prestaged )
// *****************************************************************************
//  Compute internal surface flux integrals for const-order DG (not p-adaptive)
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] t Physical time
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
// //! \param[in] solidx Material index indicator
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
// //! \param[in] dt Delta time
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  ErrChk( inciter::g_inputdeck.get< tag::flux >() == inciter::ctr::FluxType::HLLC,
          "surfInt_constP currently only supports HLLC Riemann solver" );

  const auto& ale = inciter::g_inputdeck.get< tag::ale, tag::ale >();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  //// Determine if we have solids in our problem
  //bool haveSolid = inciter::haveSolid(nmat, solidx);

  //Assert( (nmat==1 ? riemannDeriv.empty() : true), "Non-empty Riemann "
  //        "derivative vector for single material compflow" );

  // Quadrature points
  auto ng = tk::NGfa(ndof);

  // arrays for quadrature points
  std::array< std::vector< real >, 2 > coordgp;
  std::vector< real > wgp;

  coordgp[0].resize( ng );
  coordgp[1].resize( ng );
  wgp.resize( ng );

  // get quadrature point weights and coordinates for triangle
  GaussQuadratureTri( ng, coordgp, wgp );

  // Allocate memory
  std::vector< tk::real > B_l(rdof, 1.0), B_r(rdof, 1.0);
  std::array< std::vector< real >, 2 > state;
  state[0].resize(ncomp+nprim);
  state[1].resize(ncomp+nprim);

  // Values formerly looked up per quadrature point inside flux() and update_rhs_fa()
  // Hoisted since they are loop invariant and device path can't reach g_inputdeck
  auto nsld = inciter::numSolids(nmat, solidx);
  auto nstate = ncomp+nprim;

  ErrChk( nmat <= inciter::HLLCMultiMatConstP::NMAT_MAX_FLUX,
          "surfInt_constP supports max "+std::to_string(inciter::HLLCMultiMatConstP::NMAT_MAX_FLUX)+" materials" );
  ErrChk( nstate <= inciter::HLLCMultiMatConstP::NSTATE_MAX_FLUX,
          "surfInt_constP state vector supports max "+std::to_string(inciter::HLLCMultiMatConstP::NSTATE_MAX_FLUX)+" states" );

  // Strides needed by the accumulators which index R and riemannDeriv flatly
  // Allows one body of each to serve both host and device without needing override
  const std::size_t r_nprop = R.nprop();
  const std::size_t rd_ncol = riemannDeriv.empty() ? 0 : riemannDeriv[0].size();

  // Caller-owned flux buffer reused across every face and quadrature pt
  Kokkos::Array< tk::real, inciter::HLLCMultiMatConstP::NFLX_MAX > flx;

  // ---------------------------------------------------------------------------
  // Device path. Faces are the unit of parallelism, and an internal face writes
  // to both of its adjacent elements, so R and riemannDeriv are accumulated
  // atomically (rhsAccum/rdAccum pick the atomic overloads for views).
  //
  // U, P and R are expected resident already: rhs() owns that round trip
  // (prestaged). Mesh and face data are uploaded here behind residency checks.
  // ---------------------------------------------------------------------------
  ErrChk( !ale, "surfInt_constP device path does not support ALE: "
          "evaluateMeshVelTri has no device overload and W has no device view" );
  ErrChk( dev != nullptr, "surfInt_constP requires device views" );
  auto& dv = *dev;

  auto bparam =
    inciter::g_inputdeck.get< tag::multimat, tag::intsharp_param >();

  const std::size_t m_nprop = U.nprop();
  const std::size_t p_nprop = P.nprop();
  const std::size_t geo_nprop = geoElem.nprop();
  const std::size_t gf_nprop = geoFace.nprop();
  const std::size_t nelem = fd.Esuel().size()/4;
  const std::size_t nbfac = fd.Nbfac();
  const std::size_t nface = esuf.size()/2;
  const std::size_t rd_nrow = riemannDeriv.size();

  auto exec = Kokkos::DefaultExecutionSpace();

  // Mesh residency: re-upload only when the buffer was (re)allocated or the
  // calling partition changed. Same contract as the volume kernels.
  const bool mesh_ok = meshResident( dv, inpoel, coord, geoElem, nelem, nmat );

  auto solidx_h = changeToView( solidx.data(), solidx.size() );
  if (ensureDeviceCapacity( dv.solidx, "surf_solidx_d_view", solidx.size() )
      || !mesh_ok)
    Kokkos::deep_copy( dv.solidx, solidx_h );

  auto inpoel_h = changeToView( inpoel.data(), inpoel.size() );
  if (ensureDeviceCapacity( dv.inpoel, "surf_inpoel_d_view", inpoel.size() )
      || !mesh_ok)
    Kokkos::deep_copy( dv.inpoel, inpoel_h );

  auto cx_h = changeToView( cx.data(), cx.size() );
  if (ensureDeviceCapacity( dv.cx, "surf_cx_d_view", cx.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cx, cx_h );
  auto cy_h = changeToView( cy.data(), cy.size() );
  if (ensureDeviceCapacity( dv.cy, "surf_cy_d_view", cy.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cy, cy_h );
  auto cz_h = changeToView( cz.data(), cz.size() );
  if (ensureDeviceCapacity( dv.cz, "surf_cz_d_view", cz.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cz, cz_h );

  const std::size_t geoElem_size = geoElem.getSize();
  auto geoElem_h = changeToView( geoElem.getPointer(), geoElem_size );
  if (ensureDeviceCapacity( dv.geoElem, "surf_geoElem_d_view", geoElem_size )
      || !mesh_ok)
    Kokkos::deep_copy( dv.geoElem, geoElem_h );

  // Face residency
  const bool face_ok = faceResident( dv, esuf, inpofa, geoFace );

  auto esuf_h = changeToView( esuf.data(), esuf.size() );
  if (ensureDeviceCapacity( dv.esuf, "surf_esuf_d_view", esuf.size() )
      || !face_ok)
    Kokkos::deep_copy( dv.esuf, esuf_h );

  auto inpofa_h = changeToView( inpofa.data(), inpofa.size() );
  if (ensureDeviceCapacity( dv.inpofa, "surf_inpofa_d_view", inpofa.size() )
      || !face_ok)
    Kokkos::deep_copy( dv.inpofa, inpofa_h );

  const std::size_t geoFace_size = geoFace.getSize();
  auto geoFace_h = changeToView( geoFace.getPointer(), geoFace_size );
  if (ensureDeviceCapacity( dv.geoFace, "surf_geoFace_d_view", geoFace_size )
      || !face_ok)
    Kokkos::deep_copy( dv.geoFace, geoFace_h );

  // EOS block: constant for the run, upload once
  auto eos_h = changeToView( mat_blk.data(), nmat );
  if (ensureDeviceCapacity( dv.eos, "surf_eos_d_view", nmat ))
    Kokkos::deep_copy( dv.eos, eos_h );

  // riemannDeriv arrives holding the boundary-face contributions computed on
  // the host. Flatten and upload so the kernel adds internal faces on top.
  const std::size_t rd_size = rd_nrow*rd_ncol;
  std::vector< real > rd_flat( rd_size );
  for (std::size_t row=0; row<rd_nrow; ++row)
    for (std::size_t col=0; col<rd_ncol; ++col)
      rd_flat[row*rd_ncol + col] = riemannDeriv[row][col];
  uploadStaged( exec, dv.riemannDeriv, dv.stage_rd, rd_flat.data(), rd_size,
                "surf_rd_d_view" );

  // Quadrature in device storage, kept in the closure as for the volume
  // kernels: the access is a pure broadcast across the warp.
  Kokkos::Array< Kokkos::Array< real, NQUAD_MAX >, 3 > coordgpD{};
  Kokkos::Array< real, NQUAD_MAX > wgpD{};
  for (std::size_t igp=0; igp<ng; ++igp) {
    coordgpD[0][igp] = coordgp[0][igp];
    coordgpD[1][igp] = coordgp[1][igp];
    wgpD[igp] = wgp[igp];
  }

  // Shallow copies so the closure captures views, not dv
  auto solidx_d = dv.solidx;
  auto inpoel_d = dv.inpoel;
  auto cx_d = dv.cx;
  auto cy_d = dv.cy;
  auto cz_d = dv.cz;
  auto geoElem_d = dv.geoElem;
  auto geoFace_d = dv.geoFace;
  auto esuf_d = dv.esuf;
  auto inpofa_d = dv.inpofa;
  auto U_d = dv.U;
  auto P_d = dv.P;
  auto R_d = dv.R;
  auto rd_d = dv.riemannDeriv;
  auto eos_d = dv.eos;

  Kokkos::parallel_for( "surfInt_constP",
    Kokkos::RangePolicy< execution_space >( exec, nbfac, nface ),
    KOKKOS_LAMBDA( const std::size_t f )
  {
    const std::size_t el = static_cast< std::size_t >( esuf_d(2*f) );
    const std::size_t er = static_cast< std::size_t >( esuf_d(2*f+1) );

    // Extract the element coordinates
    Kokkos::Array< Kokkos::Array< real, 3 >, 4 > coordel_l, coordel_r;
    for (std::size_t i=0; i<4; ++i) {
      const auto nl = inpoel_d(4*el+i);
      coordel_l[i][0] = cx_d(nl);
      coordel_l[i][1] = cy_d(nl);
      coordel_l[i][2] = cz_d(nl);
      const auto nr = inpoel_d(4*er+i);
      coordel_r[i][0] = cx_d(nr);
      coordel_r[i][1] = cy_d(nr);
      coordel_r[i][2] = cz_d(nr);
    }

    // Compute the determinant of Jacobian matrix
    const auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    const auto detT_r =
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // Extract the face coordinates
    Kokkos::Array< Kokkos::Array< real, 3 >, 3 > coordfa;
    for (std::size_t i=0; i<3; ++i) {
      const auto n = inpofa_d(3*f+i);
      coordfa[i][0] = cx_d(n);
      coordfa[i][1] = cy_d(n);
      coordfa[i][2] = cz_d(n);
    }

    Kokkos::Array< real, 3 > fn{{ geoFace_d(f*gf_nprop+1),
                                  geoFace_d(f*gf_nprop+2),
                                  geoFace_d(f*gf_nprop+3) }};

    Kokkos::Array< real, NDOF_MAX > B_l, B_r;
    Kokkos::Array< Kokkos::Array< real, NSTATE_MAX >, 2 > state;
    Kokkos::Array< real, inciter::HLLCMultiMatConstP::NFLX_MAX > flx;

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordfa, coordgpD );

      Kokkos::Array< real, 3 > ref_gp_l{{
        Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l }};
      Kokkos::Array< real, 3 > ref_gp_r{{
        Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r }};

      // Compute the basis functions. The host overload receives vectors
      // pre-filled with 1.0; eval_basis writes only the entries above the
      // first, so the fill is reproduced here.
      for (std::size_t i=0; i<NDOF_MAX; ++i) { B_l[i] = 1.0; B_r[i] = 1.0; }
      eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );
      eval_basis( rdof, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2], B_r );

      const auto wt = wgpD[igp] * geoFace_d(f*gf_nprop+0);

      evalPolynomialSol( intsharp, ncomp, nprim, rdof, nmat, el, rdof,
        m_nprop, p_nprop, geo_nprop, bparam, solidx_d, inpoel_d,
        cx_d, cy_d, cz_d, geoElem_d, ref_gp_l, B_l, U_d, P_d, state[0] );
      evalPolynomialSol( intsharp, ncomp, nprim, rdof, nmat, er, rdof,
        m_nprop, p_nprop, geo_nprop, bparam, solidx_d, inpoel_d,
        cx_d, cy_d, cz_d, geoElem_d, ref_gp_r, B_r, U_d, P_d, state[1] );

      // compute flux. ALE is refused above, so the mesh velocity is zero.
      auto nflx = inciter::HLLCMultiMatConstP::flux( eos_d, fn, state, 0.0,
                    nmat, nsld, ncomp, nstate, solidx_d, flx );

      // Add the surface integration term to the rhs
      update_rhs_fa_constP( ncomp, nmat, ndof, ndof, ndof, wt, fn,
                            el, er, flx, nflx, B_l, B_r, solidx_d, nsld,
                            R_d, r_nprop, rd_d, rd_ncol );
    }
  });

  // riemannDeriv /= geoElem(e,0). Formerly a host loop in rhs(), done here so
  // riemannDeriv never leaves the device. Ordered after the face kernel by
  // running on the same execution space instance.
  Kokkos::parallel_for( "surfInt_constP_rdiv",
    Kokkos::RangePolicy< execution_space >( exec, 0, rd_ncol ),
    KOKKOS_LAMBDA( const std::size_t e )
  {
    const auto vol = geoElem_d(e*geo_nprop + 0);
    for (std::size_t row=0; row<rd_nrow; ++row)
      rd_d(row*rd_ncol + e) /= vol;
  });
}

void
surfIntFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  real t,
  const std::size_t rdof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const Fields& U,
  const Fields& P,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp )
// *****************************************************************************
//  Compute internal surface flux integrals for second order FV
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] t Physical time
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] srcFlag Whether the energy source was added
//! \param[in,out] R Right-hand side vector computed
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto& localFaceId = fd.FaceLocalId();

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  // Basis functions for all face-centroids of element e
  std::array< std::vector< tk::real >, 4 > Bf_array;
  for (std::size_t i=0; i<4; ++i) {
    Bf_array[i].resize(rdof);
    eval_basis(rdof, tk::fc_coord[i][0], tk::fc_coord[i][1], tk::fc_coord[i][2], Bf_array[i]);
  }

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // face normal
    std::array< real, 3 > fn{{geoFace(f,1), geoFace(f,2), geoFace(f,3)}};

    // face centroid
    std::array< real, 3 > gp{{geoFace(f,4), geoFace(f,5), geoFace(f,6)}};

    auto f_Lid = static_cast< std::size_t >(localFaceId[2*f]);
    auto f_Rid = static_cast< std::size_t >(localFaceId[2*f+1]);

    auto ref_gp_l = tk::fc_coord[f_Lid];
    auto ref_gp_r = tk::fc_coord[f_Rid];

    // Compute the basis functions
    const auto& B_l = Bf_array[f_Lid];
    const auto& B_r = Bf_array[f_Rid];

    std::array< std::vector< real >, 2 > state;

    state[0] = evalFVSol(mat_blk, intsharp, ncomp, nprim, rdof,
      nmat, el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P, srcFlag[el]);
    state[1] = evalFVSol(mat_blk, intsharp, ncomp, nprim, rdof,
      nmat, er, inpoel, coord, geoElem, ref_gp_r, B_r, U, P, srcFlag[er]);

    //safeReco(rdof, nmat, el, er, U, state);

    Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
            "appended boundary state vector" );
    Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
            "appended boundary state vector" );

    // evaluate prescribed velocity (if any)
    auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

    // compute flux
    auto fl = flux( mat_blk, fn, state, v, 0.0 );

    // compute non-conservative terms
    std::vector< tk::real > var_riemann(nmat+1, 0.0);
    for (std::size_t k=0; k<nmat+1; ++k) var_riemann[k] = fl[ncomp+k];

    auto ncf_l = nonConservativeIntFV(nmat, rdof, el, fn, U, P, var_riemann);
    auto ncf_r = nonConservativeIntFV(nmat, rdof, er, fn, U, P, var_riemann);

    // Add the surface integration term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(el, c) -= geoFace(f,0) * (fl[c] - ncf_l[c]);
      R(er, c) += geoFace(f,0) * (fl[c] - ncf_r[c]);
    }
  }
}

void
surfIntViscousFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  const Fields& T,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp )
// *****************************************************************************
//  Compute internal surface viscous flux integrals for second order FV
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] T Vector of temperatures at recent time step
//! \param[in] srcFlag Whether the energy source was added
//! \param[in,out] R Right-hand side vector computed
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  using inciter::velocityDofIdx;

  const auto& esuf = fd.Esuf();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > B_l(rdof), B_r(rdof);
  std::array< std::vector<tk::real>, 3 > dBdx_l, dBdx_r;
  for (std::size_t i=0; i<3; ++i) {
    dBdx_l[i].resize( rdof, 0 );
    dBdx_r[i].resize( rdof, 0 );
  }

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // Extract the element coordinates
    std::array< std::array< tk::real, 3>, 4 > coordel_l {{
      {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
      {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
      {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
      {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

    std::array< std::array< tk::real, 3>, 4 > coordel_r {{
      {{ cx[ inpoel[4*er  ] ], cy[ inpoel[4*er  ] ], cz[ inpoel[4*er  ] ] }},
      {{ cx[ inpoel[4*er+1] ], cy[ inpoel[4*er+1] ], cz[ inpoel[4*er+1] ] }},
      {{ cx[ inpoel[4*er+2] ], cy[ inpoel[4*er+2] ], cz[ inpoel[4*er+2] ] }},
      {{ cx[ inpoel[4*er+3] ], cy[ inpoel[4*er+3] ], cz[ inpoel[4*er+3] ] }} }};

    // Compute the determinant of Jacobian matrix
    auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto detT_r =
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );

    // face normal
    std::array< real, 3 > fn{{geoFace(f,1), geoFace(f,2), geoFace(f,3)}};

    // face centroid
    std::array< real, 3 > gp{{geoFace(f,4), geoFace(f,5), geoFace(f,6)}};

    // In order to determine the high-order solution from the left and right
    // elements at the surface quadrature points, the basis functions from
    // the left and right elements are needed. For this, a transformation to
    // the reference coordinates is necessary, since the basis functions are
    // defined on the reference tetrahedron only.
    // The transformation relations are shown below:
    //  xi   = Jacobian( coordel[0], gp, coordel[2], coordel[3] ) / detT
    //  eta  = Jacobian( coordel[0], coordel[2], gp, coordel[3] ) / detT
    //  zeta = Jacobian( coordel[0], coordel[2], coordel[3], gp ) / detT

    std::array< tk::real, 3> ref_gp_l{
      Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
      Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };
    std::array< tk::real, 3> ref_gp_r{
      Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
      Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
      Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r };

    //Compute the basis functions
    eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );
    eval_basis( rdof, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2], B_r );

    std::array< std::vector< real >, 2 > state;

    state[0] = evalFVSol(mat_blk, intsharp, ncomp, nprim, rdof,
      nmat, el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P, srcFlag[el]);
    state[1] = evalFVSol(mat_blk, intsharp, ncomp, nprim, rdof,
      nmat, er, inpoel, coord, geoElem, ref_gp_r, B_r, U, P, srcFlag[er]);

    Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
            "appended boundary state vector" );
    Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
            "appended boundary state vector" );

    // cell averaged state for computing the diffusive flux
    std::array< std::vector< real >, 2 > cellAvgState;
    std::vector< tk::real > Bcc(rdof, 0.0);
    Bcc[0] = 1.0;

    cellAvgState[0] = evalFVSol(mat_blk, 0, ncomp, nprim, rdof,
      nmat, el, inpoel, coord, geoElem, {{0.25, 0.25, 0.25}}, Bcc, U, P,
      srcFlag[el]);
    cellAvgState[1] = evalFVSol(mat_blk, 0, ncomp, nprim, rdof,
      nmat, er, inpoel, coord, geoElem, {{0.25, 0.25, 0.25}}, Bcc, U, P,
      srcFlag[er]);

    Assert( cellAvgState[0].size() == ncomp+nprim, "Incorrect size for "
            "appended cell-averaged state vector" );
    Assert( cellAvgState[1].size() == ncomp+nprim, "Incorrect size for "
            "appended cell-averaged state vector" );

    std::array< std::vector< real >, 2 > cellAvgT{{
      std::vector<real>(nmat), std::vector<real>(nmat) }};
    for (std::size_t k=0; k<nmat; ++k) {
      cellAvgT[0][k] = T(el, k*rdof);
      cellAvgT[1][k] = T(er, k*rdof);
    }

    std::array< std::array< real, 3 >, 2 > centroids{{
      {{geoElem(el,1), geoElem(el,2), geoElem(el,3)}},
      {{geoElem(er,1), geoElem(er,2), geoElem(er,3)}} }};

    // Numerical viscous flux
    // -------------------------------------------------------------------------

    // 1. Get spatial gradient from Dubiner dofs
    auto jacInv_l =
      tk::inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    tk::eval_dBdx_p1( rdof, jacInv_l, dBdx_l );
    auto jacInv_r =
      tk::inverseJacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );
    tk::eval_dBdx_p1( rdof, jacInv_r, dBdx_r );

    // 2. Average du_i/dx_j
    std::array< std::array< real, 3 >, 3 > dudx;
    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j)
        dudx[i][j] = 0.5 * (
            dBdx_l[j][1] * P(el, velocityDofIdx(nmat,i,rdof,1))
          + dBdx_l[j][2] * P(el, velocityDofIdx(nmat,i,rdof,2))
          + dBdx_l[j][3] * P(el, velocityDofIdx(nmat,i,rdof,3))
          + dBdx_r[j][1] * P(er, velocityDofIdx(nmat,i,rdof,1))
          + dBdx_r[j][2] * P(er, velocityDofIdx(nmat,i,rdof,2))
          + dBdx_r[j][3] * P(er, velocityDofIdx(nmat,i,rdof,3)) );

    // 3. average dT/dx_j
    std::vector< std::array< real, 3 > > dTdx(nmat,
      std::array< real, 3 >{{0, 0, 0}});
    for (std::size_t k=0; k<nmat; ++k) {
      auto mark = k*rdof;
      for (std::size_t j=0; j<3; ++j) {
        dTdx[k][j] = 0.5 * (
            dBdx_l[j][1] * T(el, mark+1)
          + dBdx_l[j][2] * T(el, mark+2)
          + dBdx_l[j][3] * T(el, mark+3)
          + dBdx_r[j][1] * T(er, mark+1)
          + dBdx_r[j][2] * T(er, mark+2)
          + dBdx_r[j][3] * T(er, mark+3) );
      }
    }

    // 4. Compute flux
    auto fl = modifiedGradientViscousFlux(nmat, ncomp, fn, centroids, state,
      cellAvgState, cellAvgT, dudx, dTdx);

    // -------------------------------------------------------------------------

    // Add the surface integration term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(el, c) += geoFace(f,0) * fl[c];
      R(er, c) -= geoFace(f,0) * fl[c];
    }
  }
}

void
surfIntViscousMultiSpecies(
  std::size_t nspec,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  Fields& R )
// *****************************************************************************
//  Compute internal surface viscous flux integrals for multispecies flow
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  if (ndof == 1) {
    MultiSpeciesViscousTermsP0P1 viscousRhs( nspec, rdof );
    viscousInternalFaceInt( viscousRhs, mat_blk, ndof, inpoel, coord, fd,
      geoFace, geoElem, U, P, R );
  }

  else if (ndof == 4) {
    MultiSpeciesViscousTermsDGP1 viscousRhs( nspec, rdof );
    viscousInternalFaceIntDG( viscousRhs, mat_blk, ndof, inpoel, coord, fd,
      geoFace, geoElem, U, P, R ); // No-op
  }
  
  else
    Throw( "Viscous operators only implemented for scheme = 'p0p1' and 'dgp1'." );
}

std::vector< real >
modifiedGradientViscousFlux(
  std::size_t nmat,
  std::size_t ncomp,
  const std::array< tk::real, 3 >& fn,
  const std::array< std::array< tk::real, 3 >, 2 >& centroids,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::array< std::vector< tk::real >, 2 >& cellAvgState,
  const std::array< std::vector< tk::real >, 2 >& cellAvgT,
  const std::array< std::array< real, 3 >, 3 > dudx,
  const std::vector< std::array< real, 3 > >& dTdx )
// *****************************************************************************
//  Compute the viscous fluxes from the left and right states
//! \param[in] nmat Number of materials
//! \param[in] ncomp Number of component equations in the PDE system
//! \param[in] fn Face/Surface normal
//! \param[in] centroids Left and right cell centroids
//! \param[in] state Left and right unknown/state vector
//! \param[in] cellAvgState Left and right cell-averaged unknown/state vector
//! \param[in] cellAvgT Left and right cell-averaged temperature vector
//! \param[in] dudx Average velocity gradient tensor
//! \param[in] dTdx Average temperature gradient tensor
//! \return Numerical viscous flux using the Modified Gradient approach.
//! \details The average gradient is modified according to Weiss et al. to
//!   obtain a stable discretization (average results in unstable central
//!   central difference).
//!   Ref: Weiss, J. M., Maruszewski, J. P., & Smith, W. A. (1999). Implicit
//!   solution of preconditioned Navier-Stokes equations using algebraic
//!   multigrid. AIAA journal, 37(1), 29-36.
// *****************************************************************************
{
  using inciter::velocityDofIdx;
  using inciter::volfracDofIdx;
  using inciter::momentumIdx;
  using inciter::velocityIdx;
  using inciter::volfracIdx;
  using inciter::energyIdx;

  std::vector< real > fl(ncomp, 0);

  // 1. Modify the average gradient
  std::array< real, 3 > r_f{{
    centroids[1][0] - centroids[0][0],
    centroids[1][1] - centroids[0][1],
    centroids[1][2] - centroids[0][2] }};
  real r_mag = std::sqrt(tk::dot(r_f, r_f));
  for (std::size_t i=0; i<3; ++i)
    r_f[i] /= r_mag;

  // velocity gradient
  auto dudx_m = dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx_m[i][j] -= ( tk::dot(dudx_m[i], r_f) -
        (cellAvgState[1][ncomp+velocityIdx(nmat,i)]
        - cellAvgState[0][ncomp+velocityIdx(nmat,i)])/r_mag ) * r_f[j];
  // temperature gradient
  auto dTdx_m = dTdx;
  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t j=0; j<3; ++j)
      dTdx_m[k][j] -= ( tk::dot(dTdx_m[k], r_f) -
        (cellAvgT[1][k] - cellAvgT[0][k])/r_mag ) * r_f[j];
  }

  // 2. Compute viscous stress tensor
  std::array< real, 6 > tau;
  real mu(0.0);
  std::vector< real > alLR(nmat, 0), conduct_mat(nmat, 0);
  for (std::size_t k=0; k<nmat; ++k)
    alLR[k] = 0.5*( state[0][volfracIdx(nmat,k)] + state[1][volfracIdx(nmat,k)] );
  for (std::size_t k=0; k<nmat; ++k) {
    mu += alLR[k] * inciter::getmatprop< tag::mu >(k);
    conduct_mat[k] = inciter::getmatprop< tag::mu >(k) *
      inciter::getmatprop< tag::cv >(k) * inciter::getmatprop< tag::gamma >(k)
      / 0.71;
  }

  tau[0] = mu * ( 4.0 * dudx_m[0][0] - 2.0*(dudx_m[1][1] + dudx_m[2][2]) ) / 3.0;
  tau[1] = mu * ( 4.0 * dudx_m[1][1] - 2.0*(dudx_m[0][0] + dudx_m[2][2]) ) / 3.0;
  tau[2] = mu * ( 4.0 * dudx_m[2][2] - 2.0*(dudx_m[0][0] + dudx_m[1][1]) ) / 3.0;
  tau[3] = mu * ( dudx_m[0][1] + dudx_m[1][0] );
  tau[4] = mu * ( dudx_m[0][2] + dudx_m[2][0] );
  tau[5] = mu * ( dudx_m[1][2] + dudx_m[2][1] );

  // 3. Compute viscous flux terms
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      fl[momentumIdx(nmat, i)] += tau[inciter::stressCmp[i][j]] * fn[j];

  std::vector< std::array< real, 3 > > energyFlux(nmat, {{0.0, 0.0, 0.0}});
  std::array< real, 3 > uAvg{{
    0.5*(state[0][ncomp+velocityIdx(nmat,0)] + state[1][ncomp+velocityIdx(nmat,0)]),
    0.5*(state[0][ncomp+velocityIdx(nmat,1)] + state[1][ncomp+velocityIdx(nmat,1)]),
    0.5*(state[0][ncomp+velocityIdx(nmat,2)] + state[1][ncomp+velocityIdx(nmat,2)])
    }};

  for (std::size_t k=0; k<nmat; ++k) {
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t i=0; i<3; ++i)
        energyFlux[k][j] += uAvg[i] * tau[inciter::stressCmp[i][j]] +
          conduct_mat[k] * dTdx_m[k][j];
  }

  for (std::size_t k=0; k<nmat; ++k) {
    fl[energyIdx(nmat, k)] = alLR[k]
      * tk::dot(energyFlux[k], fn);
  }

  return fl;
}

} // tk::
