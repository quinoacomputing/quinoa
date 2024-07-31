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
#include "Vector.hpp"
#include "Quadrature.hpp"
#include "Reconstruction.hpp"
#include "Integrate/SolidTerms.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "MultiMat/MiscMultiMatFns.hpp"
#include "EoS/GetMatProp.hpp"

namespace inciter {
extern ctr::InputDeck g_inputdeck;
}

namespace tk {

void
surfInt( std::size_t nmat,
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
         const std::vector< std::size_t >& ndofel,
         const tk::real /*dt*/,
         Fields& R,
         std::vector< std::vector< tk::real > >&,
         std::vector< std::vector< tk::real > >&,
         std::vector< std::vector< tk::real > >& riemannDeriv,
         int intsharp )
// *****************************************************************************
//  Compute internal surface flux integrals
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
//! \param[in] ndofel Vector of local number of degrees of freedom
// //! \param[in] dt Delta time
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] vriem Vector of the riemann velocity
//! \param[in,out] riemannLoc Vector of coordinates where Riemann velocity data
//!   is available
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
//! \param[in] intsharp Interface compression tag, an optional argument, with
//!   default 0, so that it is unused for single-material and transport.
// *****************************************************************************
{
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

      std::array< tk::real, 3> ref_gp_l{
        Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
        Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l };
      std::array< tk::real, 3> ref_gp_r{
        Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
        Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r };

      //Compute the basis functions
      auto B_l = eval_basis( dof_el, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2] );
      auto B_r = eval_basis( dof_er, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2] );

      auto wt = wgp[igp] * geoFace(f,0);

      std::array< std::vector< real >, 2 > state;

      state[0] = evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof,
        nmat, el, dof_el, inpoel, coord, geoElem, ref_gp_l, B_l, U, P);
      state[1] = evalPolynomialSol(mat_blk, intsharp, ncomp, nprim, rdof,
        nmat, er, dof_er, inpoel, coord, geoElem, ref_gp_r, B_r, U, P);

      Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );
      Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );

      // evaluate prescribed velocity (if any)
      auto v = vel( ncomp, gp[0], gp[1], gp[2], t );

      // compute flux
      auto fl = flux( mat_blk, fn, state, v );

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
    for(std::size_t idof = 0; idof < ndof; idof++) {
      riemannDeriv[3*nmat+idof][el] += wt * fl[ncomp+nmat] * B_l[idof];
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
  }
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

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

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
    auto B_l = eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2] );
    auto B_r = eval_basis( rdof, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2] );

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
    auto fl = flux( mat_blk, fn, state, v );

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
    auto B_l = eval_basis( rdof, ref_gp_l[0], ref_gp_l[1], ref_gp_l[2] );
    auto B_r = eval_basis( rdof, ref_gp_r[0], ref_gp_r[1], ref_gp_r[2] );

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

    std::array< std::array< real, 3 >, 2 > centroids{{
      {{geoElem(el,1), geoElem(el,2), geoElem(el,3)}},
      {{geoElem(er,1), geoElem(er,2), geoElem(er,3)}} }};

    // Numerical viscous flux
    // -------------------------------------------------------------------------

    // 1. Get spatial gradient from Dubiner dofs
    auto jacInv_l =
      tk::inverseJacobian( coordel_l[0], coordel_l[1], coordel_l[2], coordel_l[3] );
    auto dBdx_l = tk::eval_dBdx_p1( rdof, jacInv_l );
    auto jacInv_r =
      tk::inverseJacobian( coordel_r[0], coordel_r[1], coordel_r[2], coordel_r[3] );
    auto dBdx_r = tk::eval_dBdx_p1( rdof, jacInv_r );

    std::array< std::array< real, 3 >, 3 > dudx;

    // 2. Average du_i/dx_j
    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j)
        dudx[i][j] = 0.5 * (
            dBdx_l[j][1] * P(el, velocityDofIdx(nmat,i,rdof,1))
          + dBdx_l[j][2] * P(el, velocityDofIdx(nmat,i,rdof,2))
          + dBdx_l[j][3] * P(el, velocityDofIdx(nmat,i,rdof,3))
          + dBdx_r[j][1] * P(er, velocityDofIdx(nmat,i,rdof,1))
          + dBdx_r[j][2] * P(er, velocityDofIdx(nmat,i,rdof,2))
          + dBdx_r[j][3] * P(er, velocityDofIdx(nmat,i,rdof,3)) );

    // 3. Compute flux
    auto fl = modifiedGradientViscousFlux(nmat, ncomp, fn, centroids, state,
      cellAvgState, dudx);

    // -------------------------------------------------------------------------

    // Add the surface integration term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(el, c) += geoFace(f,0) * fl[c];
      R(er, c) -= geoFace(f,0) * fl[c];
    }
  }
}

std::vector< real >
modifiedGradientViscousFlux(
  std::size_t nmat,
  std::size_t ncomp,
  const std::array< tk::real, 3 >& fn,
  const std::array< std::array< tk::real, 3 >, 2 >& centroids,
  const std::array< std::vector< tk::real >, 2 >& state,
  const std::array< std::vector< tk::real >, 2 >& cellAvgState,
  const std::array< std::array< real, 3 >, 3 > dudx )
// *****************************************************************************
//  Compute the viscous fluxes from the left and right states
//! \param[in] nmat Number of materials
//! \param[in] ncomp Number of component equations in the PDE system
//! \param[in] fn Face/Surface normal
//! \param[in] centroids Left and right cell centroids
//! \param[in] state Left and right unknown/state vector
//! \param[in] cellAvgState Left and right cell-averaged unknown/state vector
//! \param[in] dudx Average velocity gradient tensor
//! \return Numerical viscous flux using the Modified Gradient approach.
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

  auto dudx_m = dudx;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      dudx_m[i][j] -= ( tk::dot(dudx_m[i], r_f) -
        (cellAvgState[1][ncomp+velocityIdx(nmat,i)]
        - cellAvgState[0][ncomp+velocityIdx(nmat,i)])/r_mag ) * r_f[j];

  // 2. Compute viscous stress tensor
  std::array< real, 6 > tau;
  real mu(0.0);
  std::vector< real > alLR(nmat, 0);
  for (std::size_t k=0; k<nmat; ++k)
    alLR[k] = 0.5*( state[0][volfracIdx(nmat,k)] + state[1][volfracIdx(nmat,k)] );
  for (std::size_t k=0; k<nmat; ++k)
    mu += alLR[k] * inciter::getmatprop< tag::mu >(k);

  tau[0] = mu * ( 4.0 * dudx_m[0][0] - 2.0*(dudx_m[1][1] + dudx_m[2][2]) ) / 3.0;
  tau[1] = mu * ( 4.0 * dudx_m[0][0] - 2.0*(dudx_m[1][1] + dudx_m[2][2]) ) / 3.0;
  tau[2] = mu * ( 4.0 * dudx_m[0][0] - 2.0*(dudx_m[1][1] + dudx_m[2][2]) ) / 3.0;
  tau[3] = mu * ( dudx_m[0][1] + dudx_m[1][0] );
  tau[4] = mu * ( dudx_m[0][2] + dudx_m[2][0] );
  tau[5] = mu * ( dudx_m[1][2] + dudx_m[2][1] );

  // 3. Compute viscous flux terms
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      fl[momentumIdx(nmat, i)] += tau[inciter::stressCmp[i][j]] * fn[j];

  std::array< real, 3 > energyFlux{{0.0, 0.0, 0.0}},
    uAvg{{
    0.5*(state[0][ncomp+velocityIdx(nmat,0)] + state[1][ncomp+velocityIdx(nmat,0)]),
    0.5*(state[0][ncomp+velocityIdx(nmat,1)] + state[1][ncomp+velocityIdx(nmat,1)]),
    0.5*(state[0][ncomp+velocityIdx(nmat,2)] + state[1][ncomp+velocityIdx(nmat,2)])
    }};
  for (std::size_t j=0; j<3; ++j)
    for (std::size_t i=0; i<3; ++i)
      energyFlux[j] += uAvg[i] * tau[inciter::stressCmp[i][j]];

  for (std::size_t k=0; k<nmat; ++k) {
    fl[energyIdx(nmat, k)] = alLR[k]
      * tk::dot(energyFlux, fn);
  }

  return fl;
}

} // tk::
