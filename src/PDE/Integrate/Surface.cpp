// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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

void
tk::surfInt( ncomp_t system,
             std::size_t nmat,
             ncomp_t offset,
             const std::size_t ndof,
             const std::size_t rdof,
             const std::vector< std::size_t >& inpoel,
             const UnsMesh::Coords& coord,
             const inciter::FaceData& fd,
             const Fields& geoFace,
             const RiemannFluxFn& flux,
             const VelFn& vel,
             const Fields& U,
             const Fields& P,
             const std::vector< std::size_t >& ndofel,
             Fields& R,
             std::vector< std::vector< tk::real > >& riemannDeriv )
// *****************************************************************************
//  Compute internal surface flux integrals
//! \param[in] system Equation system index
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] P Vector of primitives at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  auto ncomp = U.nprop()/rdof;
  auto nprim = P.nprop()/rdof;

  Assert( (nmat==1 ? riemannDeriv.empty() : true), "Non-empty Riemann "
          "derivative vector for single material compflow" );

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
      fn{{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }};

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

      //Compute the basis functions
      auto B_l = eval_basis( dof_el,
            Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l );
      auto B_r = eval_basis( dof_er,
            Jacobian( coordel_r[0], gp, coordel_r[2], coordel_r[3] ) / detT_r,
            Jacobian( coordel_r[0], coordel_r[1], gp, coordel_r[3] ) / detT_r,
            Jacobian( coordel_r[0], coordel_r[1], coordel_r[2], gp ) / detT_r );

      auto wt = wgp[igp] * geoFace(f,0,0);

      std::array< std::vector< real >, 2 > state;
      std::array< std::vector< real >, 2 > sprim;

      state[0] = eval_state( ncomp, offset, rdof, dof_el, el, U, B_l );
      sprim[0] = eval_state( nprim, offset, rdof, dof_el, el, P, B_l );
      state[1] = eval_state( ncomp, offset, rdof, dof_er, er, U, B_r );
      sprim[1] = eval_state( nprim, offset, rdof, dof_er, er, P, B_r );

      // consolidate primitives into state vector
      state[0].insert(state[0].end(), sprim[0].begin(), sprim[0].end());
      state[1].insert(state[1].end(), sprim[1].begin(), sprim[1].end());

      Assert( state[0].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );
      Assert( state[1].size() == ncomp+nprim, "Incorrect size for "
              "appended boundary state vector" );

      // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );

      // compute flux
      auto fl =
         flux( fn, state, v );

      // Add the surface integration term to the rhs
      update_rhs_fa( ncomp, nmat, offset, ndof, ndofel[el], ndofel[er], wt, fn,
                     el, er, fl, B_l, B_r, R, riemannDeriv );
    }
  }
}

void
tk::update_rhs_fa ( ncomp_t ncomp,
                    std::size_t nmat,
                    ncomp_t offset,
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
//! \param[in] offset Offset this PDE system operates from
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

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(el, mark, offset) -= wt * fl[c];
    R(er, mark, offset) += wt * fl[c];

    if(ndof_l > 1)          //DG(P1)
    {
      R(el, mark+1, offset) -= wt * fl[c] * B_l[1];
      R(el, mark+2, offset) -= wt * fl[c] * B_l[2];
      R(el, mark+3, offset) -= wt * fl[c] * B_l[3];
    }

    if(ndof_r > 1)          //DG(P1)
    {
      R(er, mark+1, offset) += wt * fl[c] * B_r[1];
      R(er, mark+2, offset) += wt * fl[c] * B_r[2];
      R(er, mark+3, offset) += wt * fl[c] * B_r[3];
    }

    if(ndof_l > 4)          //DG(P2)
    {
      R(el, mark+4, offset) -= wt * fl[c] * B_l[4];
      R(el, mark+5, offset) -= wt * fl[c] * B_l[5];
      R(el, mark+6, offset) -= wt * fl[c] * B_l[6];
      R(el, mark+7, offset) -= wt * fl[c] * B_l[7];
      R(el, mark+8, offset) -= wt * fl[c] * B_l[8];
      R(el, mark+9, offset) -= wt * fl[c] * B_l[9];
    }

    if(ndof_r > 4)          //DG(P2)
    {
      R(er, mark+4, offset) += wt * fl[c] * B_r[4];
      R(er, mark+5, offset) += wt * fl[c] * B_r[5];
      R(er, mark+6, offset) += wt * fl[c] * B_r[6];
      R(er, mark+7, offset) += wt * fl[c] * B_r[7];
      R(er, mark+8, offset) += wt * fl[c] * B_r[8];
      R(er, mark+9, offset) += wt * fl[c] * B_r[9];
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

    // Divergence of velocity
    riemannDeriv[3*nmat][el] += wt * fl[ncomp+nmat];
    riemannDeriv[3*nmat][er] -= wt * fl[ncomp+nmat];
  }
}
