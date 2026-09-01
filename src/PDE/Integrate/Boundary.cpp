// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Boundary.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing physical boundary surface integrals of a
     system of PDEs in DG methods
  \details   This file contains functionality for computing physical boundary
     surface integrals of a system of PDEs used in discontinuous Galerkin
     methods for various orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Basis.hpp"
#include "Boundary.hpp"
#include "ViscousTerms.hpp"
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "MultiMatTerms.hpp"
#include "MultiMat/MultiMatIndexing.hpp"
#include "Reconstruction.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "RhsAccum.hpp"
#include "MultiMat/BCFunctionsDev.hpp"
#include "Riemann/HLLCMultiMatConstP.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace {

using UnManagedMem = Kokkos::MemoryTraits<Kokkos::Unmanaged>;

// Wrap a host raw pointer in an unmanaged host view. Same helper as in
// Surface.cpp and MultiMatTerms.cpp, which keep it in anonymous namespaces.
template< typename T >
auto changeToView( T* object, std::size_t n ) {
  Kokkos::View< T*, Kokkos::LayoutLeft, Kokkos::HostSpace, UnManagedMem >
    object_view( object, n );
  return object_view;
}

} // anonymous namespace

namespace tk {

namespace {

template< class ViscousTerms >
void
viscousBoundaryFaceInt(
  const ViscousTerms& viscousRhs,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::vector< std::size_t >& bcconfig,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  Fields& R )
// *****************************************************************************
//! \brief Shared PDE-nonspecific boundary traversal for viscous surface operators
//! \tparam ViscousTerms Policy type that computes PDE-specific viscous RHS
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Number of active solution degrees of freedom
//! \param[in] bcconfig Boundary configuration vector for multiple side sets
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] t Physical time
//! \param[in] state Boundary state function
//! \param[in] gradFn Boundary gradient function
//! \param[in,out] R Right-hand side vector computed
//! \details This routine mirrors viscousInternalFaceInt() for physical
//!   boundaries. The ghost-side states are provided by the supplied StateFn.
// *****************************************************************************
{
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  const auto ncomp = static_cast< ncomp_t >( R.nprop()/ndof );

  std::array< std::vector< tk::real >, 2 > B;
  std::array< std::vector< tk::real >, 3 > dBdx_l;
  std::array< std::array< std::array< tk::real, 3 >, 4>, 2 > grad;
  std::vector< tk::real > fl( ncomp, 0.0 );

  for (const auto& s : bcconfig) {  // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        auto ndof_l = viscousRhs.localDof( el );

        // Extract the element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
          {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
          {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
          {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
          {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

        // face normal
        std::array< real, 3 > fn{{geoFace(f,1), geoFace(f,2), geoFace(f,3)}};

        // face centroid
        std::array< real, 3 > gp{{geoFace(f,4), geoFace(f,5), geoFace(f,6)}};

        // Quadrature point in element-reference-frame
        std::array< tk::real, 3> ref_gp_l{{
          Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
          Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l }};

        // Compute the basis functions for the left element
        B[0].resize(ndof_l);
        eval_basis( ndof_l, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B[0] );

        // Boundary condition, i.e. ghost state
        auto ghostState = state( ncomp, mat_blk,
          viscousRhs.stateAt(mat_blk, U, P, el, ndof_l, B[0]),
          gp[0], gp[1], gp[2], t, fn );

        // Cell centroids- [0]: left cell, [1]: ghost cell
        // The ghost-cell is a 'reflection' of the boundary cell about the
        // boundary-face. i.e. the vector pointing from the internal-cell
        // centroid to the ghost-cell centroid is normal to the face (aligned
        // with the face-normal), and has length 2*d. d is the distance between
        // the internal-cell centroid and the boundary-face. Based on this
        // information, the centroid of the ghost-cell can be computed using
        // vector algebra.
        std::array< std::array< real, 3 >, 2 > centroids;
        centroids[0] = {{geoElem(el,1), geoElem(el,2), geoElem(el,3)}};
        tk::real d = std::abs( tk::dot(fn,centroids[0]) + tk::dot(fn,gp) ) /
          std::sqrt(tk::dot(fn,fn));
        for (std::size_t i=0; i<3; ++i)
          centroids[1][i] = centroids[0][i] + 2.0*d*fn[i];

        // Boundary condition at ghost cell-center
        std::vector< tk::real > Bcc(ndof_l, 0.0);
        Bcc[0] = 1.0;
        auto ucc = viscousRhs.stateAt( mat_blk, U, P, el, ndof_l, Bcc );
        auto ghostCellAvgState = state( ncomp, mat_blk, ucc, centroids[1][0],
          centroids[1][1], centroids[1][2], t, fn );

        // Gradients of basis functions
        for (std::size_t i=0; i<3; ++i)
          dBdx_l[i].assign( ndof_l, 0.0 );
        auto jacInv_l =
          inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
        eval_dBdx_p1( ndof_l, jacInv_l, dBdx_l );

        // Compute gradients
        viscousRhs.gradientIntElem( U, P, el, dBdx_l, grad[0] );

        // Apply BCs on gradients
        std::vector< tk::real > dqdx_l(4*3,0.0);
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            dqdx_l[3*i+j] = grad[0][i][j];  // velocity gradients

        for (std::size_t j=0; j<3; ++j)
          dqdx_l[3*3+j] = grad[0][3][j];  // temperature gradients

        auto gradBC = gradFn( 4, mat_blk, dqdx_l, gp[0], gp[1], gp[2], t, fn );
        // store BC gradients into gradient vector
        for (std::size_t i=0; i<4; ++i)
          for (std::size_t j=0; j<3; ++j) {
            grad[1][i][j] = gradBC[1][3*i+j];
          }

        // Compute viscous fluxes
        viscousRhs.interiorFlux( mat_blk, ncomp, ghostState, ghostCellAvgState,
          fn, centroids, grad, fl );

        // Contribute fluxes to RHS
        for (ncomp_t c=0; c<ncomp; ++c)
          R(el, c) += geoFace(f,0) * fl[c];
      }
    }
  }
}

template< class ViscousTerms >
void
viscousBoundaryFaceIntDG(
  const ViscousTerms& viscousRhs,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::vector< std::size_t >& bcconfig,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const Fields& U,
  const Fields& P,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  Fields& R )
// *****************************************************************************
//! \brief Compute boundary surface flux integrals for a given boundary type for
//    viscous DG
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \tparam ViscousTerms Policy type that computes PDE-specific viscous RHS
//! \param[in] viscousRhs PDE-specific viscous residual policy
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Number of active solution degrees of freedom
//! \param[in] bcconfig Boundary configuration vector for multiple side sets
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] t Physical time
//! \param[in] state Boundary state function
//! \param[in] gradFn Boundary gradient function
//! \param[in,out] R Right-hand side vector computed
//! \details This routine mirrors viscousInternalFaceIntDG() for physical
//!   boundaries. The ghost-side states are provided by the supplied StateFn.
// *****************************************************************************
{
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  const auto ncomp = static_cast< ncomp_t >( R.nprop()/ndof );

  std::array < std::vector< tk::real >, 2 > B;
  std::array < std::vector< tk::real >, 3 > dBdx_l;
  std::array < std::array< std::array< tk::real, 3 >, 4>, 2 > grad;
  std::vector< tk::real > fl( ncomp, 0.0);

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);
        auto ndof_l = viscousRhs.localDof( el );
        auto ng = tk::NGfa(ndof_l);

        // arrays for quadrature points
        std::array< std::vector< real >, 2 > coordgp;
        std::vector< real > wgp;

        coordgp[0].resize( ndof_l );
        coordgp[1].resize( ndof_l );
        wgp.resize( ndof_l );

        // get quadrature point weights and coordinates for triangle
        GaussQuadratureTri( ndof_l, coordgp, wgp );

        // Extract the left element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
        {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
        {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
        {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
        {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

        // Extract the face coordinates
        std::array< std::array< tk::real, 3>, 3 > coordfa {{
          {{ cx[ inpofa[3*f  ] ], cy[ inpofa[3*f  ] ], cz[ inpofa[3*f  ] ] }},
          {{ cx[ inpofa[3*f+1] ], cy[ inpofa[3*f+1] ], cz[ inpofa[3*f+1] ] }},
          {{ cx[ inpofa[3*f+2] ], cy[ inpofa[3*f+2] ], cz[ inpofa[3*f+2] ] }} }};

        // face normal
        std::array< real, 3 > fn{{geoFace(f,1), geoFace(f,2), geoFace(f,3)}};

        // Gaussian quadrature
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the coordinates of quadrature point at physical domain
          auto gp = eval_gp( igp, coordfa, coordgp );

          std::array< tk::real, 3> ref_gp_l{
            Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };

          // Compute the basis functions for the left element
          std::vector< tk::real > B_l( ndof_l );
          eval_basis( ndof_l, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );

          auto wt = wgp[igp] * geoFace(f,0);

          // Compute the state variables at the left element
          //auto state_l = viscousRhs.stateAt( mat_blk, U, P, el, ndof_l, B_l);

          auto State = state( ncomp, mat_blk,
            viscousRhs.stateAt(mat_blk, U, P, el, ndof_l, B[0]),
            gp[0], gp[1], gp[2], t, fn );

          // Gradients of basis functions
          for (std::size_t i=0; i<3; ++i)
            dBdx_l[i].assign( ndof_l, 0.0 );
          auto jacInv_l =
            inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
          eval_dBdx_p1( ndof_l, jacInv_l, dBdx_l );

          // Compute gradients
          viscousRhs.gradientIntElem( U, P, el, dBdx_l, grad[0] ); // no-op for now

          // Apply BCs on gradients
          std::vector< tk::real > dqdx_l(4*3,0.0);
          for (std::size_t i=0; i<3; ++i)
            for (std::size_t j=0; j<3; ++j)
              dqdx_l[3*i+j] = grad[0][i][j];  // velocity gradients

          for (std::size_t j=0; j<3; ++j)
            dqdx_l[3*3+j] = grad[0][3][j];  // temperature gradients

          auto gradBC = gradFn( 4, mat_blk, dqdx_l, gp[0], gp[1], gp[2], t, fn );
          // store BC gradients into gradient vector
          for (std::size_t i=0; i<4; ++i)
            for (std::size_t j=0; j<3; ++j) {
              grad[1][i][j] = gradBC[1][3*i+j];
            }

          // Compute viscous fluxes
          viscousRhs.interiorFlux( mat_blk, ncomp, State,
            fn, grad, fl ); // no-op for now

          // Contribute fluxes to RHS
          for (ncomp_t c=0; c<ncomp; ++c)
          {
            auto mark = c*ndof;
            R(el, mark) += wt * fl[c];

            if (ndof_l > 1) // DG(P1)
            {
              R(el, mark+1) += wt * fl[c] * B_l[1];
              R(el, mark+1) += wt * fl[c] * B_l[2];
              R(el, mark+1) += wt * fl[c] * B_l[3];
            }

            if(ndof_l > 4)  //DG(P2)
            {
              R(el, mark+4) += wt * fl[c] * B_l[4];
              R(el, mark+5) += wt * fl[c] * B_l[5];
              R(el, mark+6) += wt * fl[c] * B_l[6];
              R(el, mark+7) += wt * fl[c] * B_l[7];
              R(el, mark+8) += wt * fl[c] * B_l[8];
              R(el, mark+9) += wt * fl[c] * B_l[9];
            }

          }

        }
      }
    }
  }
}

} // anonymous namespace

void
bndSurfInt( const bool pref,
            std::size_t nmat,
            const std::vector< inciter::EOS >& mat_blk,
            const std::size_t ndof,
            const std::size_t rdof,
            const std::vector< std::size_t >& bcconfig,
            const inciter::FaceData& fd,
            const Fields& geoFace,
            const Fields& geoElem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            real t,
            const RiemannFluxFn& flux,
            const VelFn& vel,
            const StateFn& state,
            const Fields& U,
            const Fields& P,
            const Fields& W,
            const std::vector< std::size_t >& ndofel,
            Fields& R,
            std::vector< std::vector< tk::real > >& riemannDeriv,
            int intsharp )
// *****************************************************************************
//! Compute boundary surface flux integrals for a given boundary type for DG
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] pref Indicator for p-adaptive algorithm
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
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
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > ugp(ncomp+nprim);

  //Assert( (nmat==1 ? riemannDeriv.empty() : true), "Non-empty Riemann "
  //        "derivative vector for single material compflow" );

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        auto ng = tk::NGfa(ndofel[el]);

        // arrays for quadrature points
        std::array< std::vector< real >, 2 > coordgp;
        std::vector< real > wgp;

        coordgp[0].resize( ng );
        coordgp[1].resize( ng );
        wgp.resize( ng );

        // get quadrature point weights and coordinates for triangle
        GaussQuadratureTri( ng, coordgp, wgp );

        // Extract the left element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
        {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
        {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
        {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
        {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

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

          // If an rDG method is set up (P0P1), then, currently we compute the P1
          // basis functions and solutions by default. This implies that P0P1 is
          // unsupported in the p-adaptive DG (PDG). This is a workaround until
          // we have rdofel, which is needed to distinguish between ndofs and
          // rdofs per element for pDG.
          std::size_t dof_el;
          if (rdof > ndof)
          {
            dof_el = rdof;
          }
          else
          {
            dof_el = ndofel[el];
          }

          // For multi-material p-adaptive simulation, when dofel = 1, p0p1 is applied and ndof
          // for solution evaluation should be 4
          if(ncomp > 5 && dof_el == 1 && pref)
            dof_el = 4;

          std::array< tk::real, 3> ref_gp_l{
            Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };

          //Compute the basis functions for the left element
          std::vector< tk::real > B_l(dof_el);
          eval_basis( dof_el, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );

          auto wt = wgp[igp] * geoFace(f,0);

          // Compute the state variables at the left element
          evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof, nmat, el,
            dof_el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P, ugp);

          auto var = state( ncomp, mat_blk, ugp, gp[0], gp[1], gp[2], t, fn );

          // mesh velocity at quadrature point
          tk::real wn_igp(0.0);
          if (ale) {
            auto w_igp = evaluateMeshVelTri( f, igp, inpofa, coordgp, W );

            // mesh velocity normal to element face
            wn_igp = tk::dot(w_igp, fn);
          }

          // Compute the numerical flux
          auto fl = flux(mat_blk, fn, var, vel(ncomp, gp[0], gp[1], gp[2], t),
            wn_igp);

          // Code below commented until details about the form of these terms in
          // the \alpha_k g_k equations are sorted out.
          // // Add RHS inverse deformation terms if necessary
          // if (haveSolid)
          //   solidTermsSurfInt( nmat, ndof, rdof, fn, el, er, solidx, geoElem, U,
          //                      coordel_l, coordel_r, igp, coordgp, dt, fl );

          // Add the surface integration term to the rhs
          update_rhs_bc( ncomp, nmat, ndof, ndofel[el], wt, fn, el, fl,
                         B_l, R, riemannDeriv );
        }
      }
    }
  }
}

void
update_rhs_bc ( ncomp_t ncomp,
                std::size_t nmat,
                const std::size_t ndof,
                const std::size_t ndof_l,
                const tk::real wt,
                const std::array< tk::real, 3 >& fn,
                const std::size_t el,
                const std::vector< tk::real >& fl,
                const std::vector< tk::real >& B_l,
                Fields& R,
                std::vector< std::vector< tk::real > >& riemannDeriv )
// *****************************************************************************
//  Update the rhs by adding the boundary surface integration term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_l Number of degrees of freedom for the left element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] fn Face/Surface normal
//! \param[in] el Left element index
//! \param[in] fl Surface flux
//! \param[in] B_l Basis function for the left element
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  // following line commented until rdofel is made available.
  //Assert( B_l.size() == ndof_l, "Size mismatch" );

  using inciter::newSolidsAccFn;

  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(el, mark) -= wt * fl[c];

    if(ndof_l > 1)          //DG(P1)
    {
      R(el, mark+1) -= wt * fl[c] * B_l[1];
      R(el, mark+2) -= wt * fl[c] * B_l[2];
      R(el, mark+3) -= wt * fl[c] * B_l[3];
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
  }

  // Prep for non-conservative terms in multimat
  if (fl.size() > ncomp)
  {
    // Gradients of partial pressures
    for (std::size_t k=0; k<nmat; ++k)
    {
      for (std::size_t idir=0; idir<3; ++idir)
        riemannDeriv[3*k+idir][el] += wt * fl[ncomp+k] * fn[idir];
    }

    // Divergence of velocity multiples basis fucntion( d(uB) / dx )
    for(std::size_t idof = 0; idof < ndof_l; idof++)
      riemannDeriv[3*nmat+idof][el] += wt * fl[ncomp+nmat] * B_l[idof];

    // Divergence of asigma: d(asig_ij)/dx_j
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
      {
        std::size_t mark = ncomp+nmat+1+3*(solidx[k]-1);

        for (std::size_t i=0; i<3; ++i)
          riemannDeriv[3*nmat+ndof+3*(solidx[k]-1)+i][el] -=
            wt * fl[mark+i];
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
              }
  }
}

//! \details Copy of update_rhs_bc, templated on its container types so one body
//!   serves the host loop and the device kernel. solidx and nsld are passed in
//!   rather than read from the input deck (a lookup per quadrature point, and
//!   unreachable from device code), and the logical length of fl arrives as
//!   nflx, since a fixed-size buffer's size() is its capacity. update_rhs_bc
//!   itself is untouched, so the p-adaptive path is unaffected.
//! \note A boundary face writes to one element only, so no atomics are
//!   required here; the accumulators still route through rhsAccum/rdAccum so
//!   that one body serves host containers and device views alike.
template< class FnT, class FlT, class BT, class SolidxT, class RT,
          class RiemannDerivT >
KOKKOS_INLINE_FUNCTION void
update_rhs_bc_constP ( ncomp_t ncomp,
                std::size_t nmat,
                const std::size_t ndof,
                const std::size_t ndof_l,
                const tk::real wt,
                const FnT& fn,
                const std::size_t el,
                const FlT& fl,
                const std::size_t nflx,
                const BT& B_l,
                const SolidxT& solidx,
                const std::size_t nsld,
                RT& R,
                const std::size_t r_nprop,
                RiemannDerivT& riemannDeriv,
                const std::size_t rd_ncol )
// *****************************************************************************
//  Update the rhs by adding the boundary surface integration term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] ndof_l Number of degrees of freedom for the left element
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] fn Face/Surface normal
//! \param[in] el Left element index
//! \param[in] fl Surface flux
//! \param[in] B_l Basis function for the left element
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  // following line commented until rdofel is made available.
  //Assert( B_l.size() == ndof_l, "Size mismatch" );

  using inciter::newSolidsAccFn;


  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    rhsAccum( R, el, r_nprop, mark, -wt * fl[c] );

    if(ndof_l > 1)          //DG(P1)
    {
      rhsAccum( R, el, r_nprop, mark+1, -wt * fl[c] * B_l[1] );
      rhsAccum( R, el, r_nprop, mark+2, -wt * fl[c] * B_l[2] );
      rhsAccum( R, el, r_nprop, mark+3, -wt * fl[c] * B_l[3] );
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
  }

  // Prep for non-conservative terms in multimat
  if (nflx > ncomp)
  {
    // Gradients of partial pressures
    for (std::size_t k=0; k<nmat; ++k)
    {
      for (std::size_t idir=0; idir<3; ++idir)
        rdAccum( riemannDeriv, 3*k+idir, el, rd_ncol,
            wt * fl[ncomp+k] * fn[idir] );
    }

    // Divergence of velocity multiples basis fucntion( d(uB) / dx )
    for(std::size_t idof = 0; idof < ndof_l; idof++)
      rdAccum( riemannDeriv, 3*nmat+idof, el, rd_ncol,
            wt * fl[ncomp+nmat] * B_l[idof] );

    // Divergence of asigma: d(asig_ij)/dx_j
    for (std::size_t k=0; k<nmat; ++k)
      if (solidx[k] > 0)
      {
        std::size_t mark = ncomp+nmat+1+3*(solidx[k]-1);

        for (std::size_t i=0; i<3; ++i)
          rdAccum( riemannDeriv, 3*nmat+ndof+3*(solidx[k]-1)+i, el, rd_ncol,
            -wt * fl[mark+i] );
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
              }
  }
}

void
bndSurfInt_constP(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const StateFn& state,
  const Fields& U,
  const Fields& P,
  const Fields& W,
  Fields& R,
  std::vector< std::vector< tk::real > >& riemannDeriv,
  int intsharp )
// *****************************************************************************
//! \brief Compute boundary surface flux integrals for a given boundary type for
//!   const-order DG (not p-adaptive)
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
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
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

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
  std::vector< tk::real > B_l(rdof, 1.0), ugp(ncomp+nprim);

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        // Extract the left element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
        {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
        {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
        {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
        {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

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

          //Compute the basis functions for the left element
          eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );

          auto wt = wgp[igp] * geoFace(f,0);

          // Compute the state variables at the left element
          evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof, nmat, el,
            rdof, inpoel, coord, geoElem, ref_gp_l, B_l, U, P, ugp);

          auto var = state( ncomp, mat_blk, ugp, gp[0], gp[1], gp[2], t, fn );

          // mesh velocity at quadrature point
          tk::real wn_igp(0.0);
          if (ale) {
            auto w_igp = evaluateMeshVelTri( f, igp, inpofa, coordgp, W );

            // mesh velocity normal to element face
            wn_igp = tk::dot(w_igp, fn);
          }

          // Compute the numerical flux
          auto fl = flux(mat_blk, fn, var, vel(ncomp, gp[0], gp[1], gp[2], t),
            wn_igp);

          // Code below commented until details about the form of these terms in
          // the \alpha_k g_k equations are sorted out.
          // // Add RHS inverse deformation terms if necessary
          // if (haveSolid)
          //   solidTermsSurfInt( nmat, ndof, rdof, fn, el, er, solidx, geoElem, U,
          //                      coordel_l, coordel_r, igp, coordgp, dt, fl );

          // Add the surface integration term to the rhs
          update_rhs_bc( ncomp, nmat, ndof, ndof, wt, fn, el, fl,
                         B_l, R, riemannDeriv );
        }
      }
    }
  }
}
void
bndSurfIntMultiMat_constP(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::pair< std::vector< std::size_t >, int > >& bcsets,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const Fields& U,
  const Fields& P,
  const Fields& W,
  Fields& R,
  std::vector< std::vector< tk::real > >& riemannDeriv,
  int intsharp,
  BndSurfIntDeviceViews* dev,
  [[maybe_unused]] bool prestaged )
// *****************************************************************************
//! \brief Compute boundary surface flux integrals for const-order multi-material
//!   DG (not p-adaptive)
//! \details MultiMat-specific: the Riemann solver is fixed to HLLC and the
//!   boundary state comes from inciter::bcStateDev, so the flux, vel and state
//!   function arguments the shared bndSurfInt_constP takes are not needed here.
//!   bndSurfInt_constP is left untouched for CompFlow, MultiSpecies and
//!   Transport, which pass their own state functions.
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcsets Supported boundary conditions for this mesh: for each,
//!   the side sets it applies to and its inciter::BCKind. All of them are
//!   handled in one pass so that the device version is a single launch rather
//!   than one per boundary condition.
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] W Mesh velocity vector at recent time step
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
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  // Loop-invariants for the _constP path. solidx and nsld were read from the
  // input deck inside the state function and update_rhs_bc, once per
  // quadrature point; flx replaces a per-point heap allocation.
  const auto& solidx =
    inciter::g_inputdeck.get< tag::matidxmap, tag::solidx >();
  auto nsld = inciter::numSolids(nmat, solidx);
  const std::size_t r_nprop = R.nprop();
  const std::size_t rd_ncol = riemannDeriv.empty() ? 0 : riemannDeriv[0].size();
  Kokkos::Array< tk::real, inciter::HLLCMultiMatConstP::NFLX_MAX > flx;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

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
  std::vector< tk::real > B_l(rdof, 1.0), ugp(ncomp+nprim);

  // Flatten every supported boundary condition into one face list, paired
  // with the state function each face needs. Visitation order is preserved
  // exactly: bcsets order, then side set order within each (std::map iterates
  // by id, as the previous nested loop did), then face order within each side
  // set. Preserving it keeps the accumulation into R bit-identical.
  std::vector< std::size_t > bfaces;
  std::vector< int > bfkind;
  for (const auto& bs : bcsets)
    for (const auto& s : bs.first) {
      auto bc = bface.find(static_cast<int>(s));
      if (bc != end(bface))
        for (const auto& f : bc->second) {
          bfaces.push_back( f );
          bfkind.push_back( bs.second );
        }
    }

  // ---------------------------------------------------------------------------
  // Device path. One launch covers every supported boundary condition: bfaces
  // holds the flattened face list and bfkind the state function each face
  // needs, so bcStateDev's branch is uniform within a side set.
  //
  // A boundary face writes to one element only, so no atomics are required;
  // rhsAccum/rdAccum still route through the view overloads, whose atomicity
  // is harmless here.
  //
  // U, P and R must already be resident: rhs() owns that round trip. This
  // kernel writes riemannDeriv, which rhs() zeroes on the device beforehand.
  // ---------------------------------------------------------------------------
  ErrChk( !ale, "bndSurfIntMultiMat_constP device path does not support ALE: "
          "evaluateMeshVelTri has no device overload and W has no device view" );
  ErrChk( dev != nullptr, "bndSurfIntMultiMat_constP requires device views" );
  auto& dv = *dev;

  if (bfaces.empty()) return;

  auto bparam =
    inciter::g_inputdeck.get< tag::multimat, tag::intsharp_param >();

  const std::size_t m_nprop = U.nprop();
  const std::size_t p_nprop = P.nprop();
  const std::size_t geo_nprop = geoElem.nprop();
  const std::size_t gf_nprop = geoFace.nprop();
  const std::size_t nelem = fd.Esuel().size()/4;
  const std::size_t nbf = bfaces.size();
  const std::size_t nstate = ncomp+nprim;

  auto exec = Kokkos::DefaultExecutionSpace();

  // Mesh residency
  const bool mesh_ok = meshResident( dv, inpoel, coord, geoElem, nelem, nmat );

  auto solidx_h = changeToView( solidx.data(), solidx.size() );
  if (ensureDeviceCapacity( dv.solidx, "bnd_solidx_d_view", solidx.size() )
      || !mesh_ok)
    Kokkos::deep_copy( dv.solidx, solidx_h );

  auto inpoel_h = changeToView( inpoel.data(), inpoel.size() );
  if (ensureDeviceCapacity( dv.inpoel, "bnd_inpoel_d_view", inpoel.size() )
      || !mesh_ok)
    Kokkos::deep_copy( dv.inpoel, inpoel_h );

  auto cx_h = changeToView( cx.data(), cx.size() );
  if (ensureDeviceCapacity( dv.cx, "bnd_cx_d_view", cx.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cx, cx_h );
  auto cy_h = changeToView( cy.data(), cy.size() );
  if (ensureDeviceCapacity( dv.cy, "bnd_cy_d_view", cy.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cy, cy_h );
  auto cz_h = changeToView( cz.data(), cz.size() );
  if (ensureDeviceCapacity( dv.cz, "bnd_cz_d_view", cz.size() ) || !mesh_ok)
    Kokkos::deep_copy( dv.cz, cz_h );

  const std::size_t geoElem_size = geoElem.getSize();
  auto geoElem_h = changeToView( geoElem.getPointer(), geoElem_size );
  if (ensureDeviceCapacity( dv.geoElem, "bnd_geoElem_d_view", geoElem_size )
      || !mesh_ok)
    Kokkos::deep_copy( dv.geoElem, geoElem_h );

  // Face residency
  const bool face_ok = faceResident( dv, esuf, inpofa, geoFace );

  auto esuf_h = changeToView( esuf.data(), esuf.size() );
  if (ensureDeviceCapacity( dv.esuf, "bnd_esuf_d_view", esuf.size() )
      || !face_ok)
    Kokkos::deep_copy( dv.esuf, esuf_h );

  auto inpofa_h = changeToView( inpofa.data(), inpofa.size() );
  if (ensureDeviceCapacity( dv.inpofa, "bnd_inpofa_d_view", inpofa.size() )
      || !face_ok)
    Kokkos::deep_copy( dv.inpofa, inpofa_h );

  const std::size_t geoFace_size = geoFace.getSize();
  auto geoFace_h = changeToView( geoFace.getPointer(), geoFace_size );
  if (ensureDeviceCapacity( dv.geoFace, "bnd_geoFace_d_view", geoFace_size )
      || !face_ok)
    Kokkos::deep_copy( dv.geoFace, geoFace_h );

  // EOS block: constant for the run
  const std::size_t eos_bytes = nmat*sizeof(inciter::EOS);
  auto eos_h = changeToView(
    reinterpret_cast< char* >( const_cast< inciter::EOS* >( mat_blk.data() ) ),
    eos_bytes );
  if (ensureDeviceCapacity( dv.eos, "bnd_eos_d_view", eos_bytes ))
    Kokkos::deep_copy( dv.eos, eos_h );

  // The flattened face list changes whenever the boundary configuration does,
  // so it is uploaded every call. It is O(nbfac), far smaller than the mesh.
  auto bfaces_h = changeToView( bfaces.data(), nbf );
  ensureDeviceCapacity( dv.bfaces, "bnd_bfaces_d_view", nbf );
  Kokkos::deep_copy( dv.bfaces, bfaces_h );
  auto bfkind_h = changeToView( bfkind.data(), nbf );
  ensureDeviceCapacity( dv.bfkind, "bnd_bfkind_d_view", nbf );
  Kokkos::deep_copy( dv.bfkind, bfkind_h );

  // Quadrature in device storage
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
  auto bfaces_d = dv.bfaces;
  auto bfkind_d = dv.bfkind;
  auto U_d = dv.U;
  auto P_d = dv.P;
  auto R_d = dv.R;
  auto rd_d = dv.riemannDeriv;
  auto eos_p = reinterpret_cast< const inciter::EOS* >( dv.eos.data() );

  Kokkos::parallel_for( "bndSurfIntMultiMat_constP",
    Kokkos::RangePolicy< execution_space >( exec, 0, nbf ),
    KOKKOS_LAMBDA( const std::size_t ifa )
  {
    const std::size_t f = bfaces_d(ifa);
    const auto bck = static_cast< inciter::BCKind >( bfkind_d(ifa) );
    const std::size_t el = static_cast< std::size_t >( esuf_d(2*f) );

    // Extract the left element coordinates
    Kokkos::Array< Kokkos::Array< real, 3 >, 4 > coordel_l;
    for (std::size_t i=0; i<4; ++i) {
      const auto n = inpoel_d(4*el+i);
      coordel_l[i][0] = cx_d(n);
      coordel_l[i][1] = cy_d(n);
      coordel_l[i][2] = cz_d(n);
    }

    // Compute the determinant of Jacobian matrix
    const auto detT_l =
      Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

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

    Kokkos::Array< real, NDOF_MAX > B_l;
    Kokkos::Array< Kokkos::Array< real, NSTATE_MAX >, 2 > var;
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

      // The host overload receives a vector pre-filled with 1.0; eval_basis
      // writes only the entries above the first, so reproduce the fill
      for (std::size_t i=0; i<NDOF_MAX; ++i) B_l[i] = 1.0;
      eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );

      const auto wt = wgpD[igp] * geoFace_d(f*gf_nprop+0);

      // Compute the state variables at the left element
      evalPolynomialSol( intsharp, ncomp, nprim, rdof, nmat, el, rdof,
        m_nprop, p_nprop, geo_nprop, bparam, solidx_d, inpoel_d,
        cx_d, cy_d, cz_d, geoElem_d, ref_gp_l, B_l, U_d, P_d, var[0] );

      // Ghost state
      inciter::bcStateDev( bck, ncomp, nmat, nstate, solidx_d, fn,
                           var[0], var[1] );

      // Compute the numerical flux. ALE is refused above, so wn is zero.
      auto nflx = inciter::HLLCMultiMatConstP::flux( eos_p, fn, var, 0.0,
                    nmat, nsld, ncomp, nstate, solidx_d, flx );

      // Add the surface integration term to the rhs
      update_rhs_bc_constP( ncomp, nmat, ndof, ndof, wt, fn, el, flx,
                            nflx, B_l, solidx_d, nsld,
                            R_d, r_nprop, rd_d, rd_ncol );
    }
  });
}

void
bndSurfIntFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const RiemannFluxFn& flux,
  const VelFn& vel,
  const StateFn& state,
  const Fields& U,
  const Fields& P,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp )
// *****************************************************************************
//! Compute boundary surface flux integrals for a given boundary type for FV
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] srcFlag Whether the energy source was added
//! \param[in,out] R Right-hand side vector computed
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
  const auto& bface = fd.Bface();
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

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        // face normal
        std::array< real, 3 >
          fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

        // face centroid
        std::array< real, 3 >
          gp{{ geoFace(f,4), geoFace(f,5), geoFace(f,6) }};

        auto f_Lid = static_cast< std::size_t >(localFaceId[2*f]);
        auto ref_gp_l = tk::fc_coord[f_Lid];

        // Compute the basis functions for the left element
        const auto& B_l = Bf_array[f_Lid];

        // Compute the state variables at the left element
        auto ugp = evalFVSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P,
          srcFlag[el]);

        Assert( ugp.size() == ncomp+nprim, "Incorrect size for "
                "appended boundary state vector" );

        auto var = state( ncomp, mat_blk, ugp, gp[0], gp[1], gp[2], t, fn );

        // Compute the numerical flux
        auto fl = flux( mat_blk, fn, var, vel(ncomp, gp[0], gp[1], gp[2], t),
          0.0 );

        // compute non-conservative terms
        std::vector< tk::real > var_riemann(nmat+1, 0.0);
        for (std::size_t k=0; k<nmat+1; ++k) var_riemann[k] = fl[ncomp+k];

        auto ncf_l = nonConservativeIntFV(nmat, rdof, el, fn, U, P, var_riemann);

        // Add the surface integration term to the rhs
        for (ncomp_t c=0; c<ncomp; ++c)
        {
          R(el, c) -= geoFace(f,0) * (fl[c] - ncf_l[c]);
        }
      }
    }
  }
}

void
bndSurfIntViscousFV(
  std::size_t nmat,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  const Fields& U,
  const Fields& P,
  const Fields& T,
  const std::vector< int >& srcFlag,
  Fields& R,
  int intsharp )
// *****************************************************************************
//! Compute boundary surface flux integrals for a given boundary type for FV
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] mat_blk EOS material block
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] gradFn Function to evaluate the left and right solution gradients
//!   at boundaries
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

  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  std::vector< tk::real > B_l(rdof);
  std::array< std::vector<tk::real>, 3 > dBdx_l;
  for (std::size_t i=0; i<3; ++i) dBdx_l[i].resize( rdof, 0 );

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find(static_cast<int>(s));// faces for side set
    if (bc != end(bface))
    {
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        // Extract the left element coordinates
        std::array< std::array< tk::real, 3>, 4 > coordel_l {{
        {{ cx[ inpoel[4*el  ] ], cy[ inpoel[4*el  ] ], cz[ inpoel[4*el  ] ] }},
        {{ cx[ inpoel[4*el+1] ], cy[ inpoel[4*el+1] ], cz[ inpoel[4*el+1] ] }},
        {{ cx[ inpoel[4*el+2] ], cy[ inpoel[4*el+2] ], cz[ inpoel[4*el+2] ] }},
        {{ cx[ inpoel[4*el+3] ], cy[ inpoel[4*el+3] ], cz[ inpoel[4*el+3] ] }} }};

        // Compute the determinant of Jacobian matrix
        auto detT_l =
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );

        // face normal
        std::array< real, 3 >
          fn{{ geoFace(f,1), geoFace(f,2), geoFace(f,3) }};

        // face centroid
        std::array< real, 3 >
          gp{{ geoFace(f,4), geoFace(f,5), geoFace(f,6) }};

        std::array< tk::real, 3> ref_gp_l{
          Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
          Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
          Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };

        //Compute the basis functions for the left element
        eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2], B_l );

        // Compute the state variables at the left element
        auto ugp = evalFVSol(mat_blk, intsharp, ncomp, nprim,
          rdof, nmat, el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P,
          srcFlag[el]);

        Assert( ugp.size() == ncomp+nprim, "Incorrect size for "
                "appended boundary state vector" );

        auto var = state( ncomp, mat_blk, ugp, gp[0], gp[1], gp[2], t, fn );

        // cell averaged state for computing the diffusive flux
        std::vector< tk::real > Bcc(rdof, 0.0);
        Bcc[0] = 1.0;

        auto ucc = evalFVSol(mat_blk, 0, ncomp, nprim, rdof,
          nmat, el, inpoel, coord, geoElem, {{0.25, 0.25, 0.25}}, Bcc, U, P,
          srcFlag[el]);

        Assert( ucc.size() == ncomp+nprim, "Incorrect size for "
                "appended cell-averaged state vector" );

        // Cell centroids- [0]: left cell, [1]: ghost cell
        // The ghost-cell is a 'reflection' of the boundary cell about the
        // boundary-face. i.e. the vector pointing from the internal-cell
        // centroid to the ghost-cell centroid is normal to the face (aligned
        // with the face-normal), and has length 2*d. d is the distance between
        // the internal-cell centroid and the boundary-face. Based on this
        // information, the centroid of the ghost-cell can be computed using
        // vector algebra.
        std::array< std::array< tk::real, 3 >, 2 > centroids;
        centroids[0] = {{geoElem(el,1), geoElem(el,2), geoElem(el,3)}};
        tk::real d = std::abs( tk::dot(fn,centroids[0]) + tk::dot(fn,gp) ) /
          std::sqrt(tk::dot(fn,fn));
        for (std::size_t i=0; i<3; ++i)
          centroids[1][i] = centroids[0][i] + 2.0*d*fn[i];

        // Get BC for cell-averaged state
        auto varcc = state( ncomp, mat_blk, ucc,
          centroids[1][0], centroids[1][1], centroids[1][2], t, fn );
        // Hard-coded temperature BC- only works for adiabatic 'noslipwall',
        // 'symmetry', 'extrapolate' bc
        std::array< std::vector< real >, 2 > Tcc{{
          std::vector<real>(nmat), std::vector<real>(nmat) }};
        for (std::size_t k=0; k<nmat; ++k) {
          Tcc[0][k] = T(el, k*rdof);
          Tcc[1][k] = Tcc[0][k];  // actual bc
        }

        // Numerical viscous flux
        // ---------------------------------------------------------------------

        // 1. Get spatial gradient from Dubiner dofs
        auto jacInv_l = tk::inverseJacobian( coordel_l[0], coordel_l[1],
          coordel_l[2], coordel_l[3] );
        tk::eval_dBdx_p1( rdof, jacInv_l, dBdx_l );

        std::vector< real > dudx_l(9,0.0);
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            dudx_l[3*i+j] =
                dBdx_l[j][1] * P(el, velocityDofIdx(nmat,i,rdof,1))
              + dBdx_l[j][2] * P(el, velocityDofIdx(nmat,i,rdof,2))
              + dBdx_l[j][3] * P(el, velocityDofIdx(nmat,i,rdof,3));

        // 2. Average du_i/dx_j
        auto grad = gradFn( 3, mat_blk, dudx_l, gp[0], gp[2], gp[2], t, fn );
        std::array< std::array< tk::real, 3 >, 3 > dudx;
        for (std::size_t i=0; i<3; ++i)
          for (std::size_t j=0; j<3; ++j)
            dudx[i][j] = 0.5 * (grad[0][3*i+j] + grad[1][3*i+j]);

        // 3. average dT/dx_j
        std::vector< std::array< real, 3 > > dTdx(nmat,
          std::array< real, 3 >{{0, 0, 0}});

        // 4. Compute flux
        auto fl = modifiedGradientViscousFlux(nmat, ncomp, fn, centroids, var,
          varcc, Tcc, dudx, dTdx);

        // Add the surface integration term to the rhs
        for (ncomp_t c=0; c<ncomp; ++c)
        {
          R(el, c) += geoFace(f,0) * fl[c];
        }
      }
    }
  }
}

void
bndSurfIntViscousMultiSpecies(
  std::size_t nspec,
  const std::vector< inciter::EOS >& mat_blk,
  const std::size_t ndof,
  const std::size_t rdof,
  const std::vector< std::size_t >& bcconfig,
  const inciter::FaceData& fd,
  const Fields& geoFace,
  const Fields& geoElem,
  const std::vector< std::size_t >& inpoel,
  const UnsMesh::Coords& coord,
  real t,
  const StateFn& state,
  const StateFn& gradFn,
  const Fields& U,
  const Fields& P,
  Fields& R )
// *****************************************************************************
//  Compute boundary surface viscous flux integrals for multispecies flow
//! \param[in] nspec Number of species in this PDE system
//! \param[in] mat_blk Material EOS block
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig Boundary configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] state Boundary state function
//! \param[in] gradFn Boundary gradient function
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  if (ndof == 1) {
    MultiSpeciesViscousTermsP0P1 viscousRhs( nspec, rdof );
    viscousBoundaryFaceInt( viscousRhs, mat_blk, ndof, bcconfig, inpoel,
      coord, fd, geoFace, geoElem, U, P, t, state, gradFn, R );
  }

  else if (ndof == 4) {
    MultiSpeciesViscousTermsDGP1 viscousRhs( nspec, rdof );
    viscousBoundaryFaceIntDG( viscousRhs, mat_blk, ndof, bcconfig, inpoel,
      coord, fd, geoFace, geoElem, U, P, t, state, gradFn, R ); // no-op
  }

  else
    Throw( "Viscous operators only implemented for scheme = 'p0p1' and 'dgp1'." );
}

} // tk::
