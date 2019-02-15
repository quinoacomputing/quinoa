// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing internal surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing internal surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Surface.h"
#include "Vector.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::surfInt( ncomp_t system,
             ncomp_t ncomp,
             ncomp_t offset,
             const std::vector< std::size_t >& inpoel,
             const UnsMesh::Coords& coord,
             const inciter::FaceData& fd,
             const Fields& geoFace,
             const RiemannFluxFn& flux,
             const VelFn& vel,
             const Fields& U,
             const Fields& limFunc,
             Fields& R )
// *****************************************************************************
//  Compute internal surface flux integrals
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  // Number of quadrature points for face integration
  auto ng = tk::NGfa(ndof);

  // arrays for quadrature points
  std::array< std::vector< real >, 2 > coordgp;
  std::vector< real > wgp;

  coordgp[0].resize( ng );
  coordgp[1].resize( ng );
  wgp.resize( ng );

  // get quadrature point weights and coordinates for triangle
  GaussQuadratureTri( ng, coordgp, wgp );

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // Nodal Coordinates of face
  std::array< std::array< real, 3>, 3 > coordfa;

  // Nodal Coordinates of element
  std::array< std::array< real, 3>, 4> coordel_l;
  std::array< std::array< real, 3>, 4> coordel_r;

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // Compute the determinant of Jacobian matrix
    auto detT_l = eval_det( el, cx, cy, cz, inpoel, coordel_l );
    auto detT_r = eval_det( er, cx, cy, cz, inpoel, coordel_r );

    coordfa[0][0] = cx[ inpofa[3*f]   ];
    coordfa[0][1] = cy[ inpofa[3*f]   ];
    coordfa[0][2] = cz[ inpofa[3*f]   ];

    coordfa[1][0] = cx[ inpofa[3*f+1] ];
    coordfa[1][1] = cy[ inpofa[3*f+1] ];
    coordfa[1][2] = cz[ inpofa[3*f+1] ];

    coordfa[2][0] = cx[ inpofa[3*f+2] ];
    coordfa[2][1] = cy[ inpofa[3*f+2] ];
    coordfa[2][2] = cz[ inpofa[3*f+2] ];

    // Gaussian quadrature
    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordfa, coordgp );

      // Compute the basis function
      auto B_l = eval_basis_fa( ndof, coordel_l, detT_l, gp );
      auto B_r = eval_basis_fa( ndof, coordel_r, detT_r, gp );

      auto wt = wgp[igp] * geoFace(f,0,0);

      std::array< std::vector< real >, 2 > state;

      state[0] = eval_state( ncomp, offset, ndof, el, U, limFunc, B_l );
      state[1] = eval_state( ncomp, offset, ndof, er, U, limFunc, B_r );

      Assert( state[0].size() == ncomp, "Size mismatch" );
      Assert( state[1].size() == ncomp, "Size mismatch" );

      // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );

      // compute flux
      auto fl =
         flux( {{geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0)}}, state, v );

      // Add the surface integration term to the rhs
      update_rhs_fa( ncomp, offset, ndof, wt, el, er, fl, B_l, B_r, R );
    }
  }
}

void
tk::update_rhs_fa ( ncomp_t ncomp,
                    ncomp_t offset,
                    const std::size_t ndof,
                    const tk::real wt,
                    const std::size_t el,
                    const std::size_t er,
                    const std::vector< tk::real >& fl,
                    const std::vector< tk::real >& B_l,
                    const std::vector< tk::real >& B_r,
                    Fields& R )
// *****************************************************************************
//  Update the rhs by adding the surface integration term
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Number of degree of freedom
//! \param[in] wt Weight of gauss quadrature point
//! \param[in] el Left element index
//! \param[in] er Right element index
//! \param[in] fl Surface flux
//! \param[in] B_l Basis function for the left element
//! \param[in] B_r Basis function for the right element
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  Assert( B_l.size() == ndof, "Size mismatch" );
  Assert( B_r.size() == ndof, "Size mismatch" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(el, mark, offset) -= wt * fl[c];
    R(er, mark, offset) += wt * fl[c];

    if(ndof > 1)          //DG(P1)
    {
      R(el, mark+1, offset) -= wt * fl[c] * B_l[1];
      R(el, mark+2, offset) -= wt * fl[c] * B_l[2];
      R(el, mark+3, offset) -= wt * fl[c] * B_l[3];

      R(er, mark+1, offset) += wt * fl[c] * B_r[1];
      R(er, mark+2, offset) += wt * fl[c] * B_r[2];
      R(er, mark+3, offset) += wt * fl[c] * B_r[3];
    }

    if(ndof > 4)          //DG(P2)
    {
      R(el, mark+4, offset) -= wt * fl[c] * B_l[4];
      R(el, mark+5, offset) -= wt * fl[c] * B_l[5];
      R(el, mark+6, offset) -= wt * fl[c] * B_l[6];
      R(el, mark+7, offset) -= wt * fl[c] * B_l[7];
      R(el, mark+8, offset) -= wt * fl[c] * B_l[8];
      R(el, mark+9, offset) -= wt * fl[c] * B_l[9];

      R(er, mark+4, offset) += wt * fl[c] * B_r[4];
      R(er, mark+5, offset) += wt * fl[c] * B_r[5];
      R(er, mark+6, offset) += wt * fl[c] * B_r[6];
      R(er, mark+7, offset) += wt * fl[c] * B_r[7];
      R(er, mark+8, offset) += wt * fl[c] * B_r[8];
      R(er, mark+9, offset) += wt * fl[c] * B_r[9];
    }
  }
}
