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
tk::surfIntP0( ncomp_t system,
               ncomp_t ncomp,
               ncomp_t offset,
               const inciter::FaceData& fd,
               const Fields& geoFace,
               const RiemannFluxFn& flux,
               const VelFn& vel,
               const Fields& U,
               Fields& R )
// *****************************************************************************
//  Compute internal surface flux integrals for DG(P0)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();

  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    auto ela = esuf[2*f];
    auto era = esuf[2*f+1];
    Assert( ela > -1 && era > -1, "Interior element detected as -1" );
    IGNORE(ela);
    IGNORE(era);

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);
    auto farea = geoFace(f,0,0);

    std::array< std::vector< real >, 2 >
      state{{ std::vector< real >( ncomp, 0.0 ),
              std::vector< real >( ncomp, 0.0 ) }};

    for (ncomp_t c=0; c<ncomp; ++c) {
      state[0][c] = U(el, c, offset);
      state[1][c] = U(er, c, offset);
    }

    // evaluate prescribed velocity (if any)
    auto v =
      vel( system, ncomp, geoFace(f,4,0), geoFace(f,5,0), geoFace(f,6,0) );

    auto fl =
       flux( {{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }}, state, v );

    for (ncomp_t c=0; c<ncomp; ++c) {
      R(el, c, offset) -= farea * fl[c];
      R(er, c, offset) += farea * fl[c];
    }
  }
}

void
tk::surfIntP1( ncomp_t system,
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
//  Compute internal surface flux integrals for DG(P1)
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

  // Number of integration points
  constexpr std::size_t NG = 3;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 2 > coordgp;
  std::array< real, NG > wgp;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  // get quadrature point weights and coordinates for triangle
  GaussQuadratureTri( coordgp, wgp );

  // compute internal surface flux integrals
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // nodal coordinates of the left element
    std::array< real, 3 >
      p1_l{{ cx[ inpoel[4*el] ],
             cy[ inpoel[4*el] ],
             cz[ inpoel[4*el] ] }},
      p2_l{{ cx[ inpoel[4*el+1] ],
             cy[ inpoel[4*el+1] ],
             cz[ inpoel[4*el+1] ] }},
      p3_l{{ cx[ inpoel[4*el+2] ],
             cy[ inpoel[4*el+2] ],
             cz[ inpoel[4*el+2] ] }},
      p4_l{{ cx[ inpoel[4*el+3] ],
             cy[ inpoel[4*el+3] ],
             cz[ inpoel[4*el+3] ] }};

    // nodal coordinates of the right element
    std::array< real, 3 >
      p1_r{{ cx[ inpoel[4*er] ],
             cy[ inpoel[4*er] ],
             cz[ inpoel[4*er] ] }},
      p2_r{{ cx[ inpoel[4*er+1] ],
             cy[ inpoel[4*er+1] ],
             cz[ inpoel[4*er+1] ] }},
      p3_r{{ cx[ inpoel[4*er+2] ],
             cy[ inpoel[4*er+2] ],
             cz[ inpoel[4*er+2] ] }},
      p4_r{{ cx[ inpoel[4*er+3] ],
             cy[ inpoel[4*er+3] ],
             cz[ inpoel[4*er+3] ] }};

    auto detT_l = Jacobian( p1_l, p2_l, p3_l, p4_l );
    auto detT_r = Jacobian( p1_r, p2_r, p3_r, p4_r );

    auto x1 = cx[ inpofa[3*f]   ];
    auto y1 = cy[ inpofa[3*f]   ];
    auto z1 = cz[ inpofa[3*f]   ];

    auto x2 = cx[ inpofa[3*f+1] ];
    auto y2 = cy[ inpofa[3*f+1] ];
    auto z2 = cz[ inpofa[3*f+1] ];

    auto x3 = cx[ inpofa[3*f+2] ];
    auto y3 = cy[ inpofa[3*f+2] ];
    auto z3 = cz[ inpofa[3*f+2] ];

    // Gaussian quadrature
    for (std::size_t igp=0; igp<NG; ++igp)
    {
      // Barycentric coordinates for the triangular face
      auto shp1 = 1.0 - coordgp[0][igp] - coordgp[1][igp];
      auto shp2 = coordgp[0][igp];
      auto shp3 = coordgp[1][igp];

      // transformation of the quadrature point from the 2D reference/master
      // element to physical domain, to obtain its physical (x,y,z)
      // coordinates.
      std::array< real, 3 > gp{{ x1*shp1 + x2*shp2 + x3*shp3,
                                 y1*shp1 + y2*shp2 + y3*shp3,
                                 z1*shp1 + z2*shp2 + z3*shp3 }};

      // The basis functions chosen for the DG method are the Dubiner
      // basis, which are Legendre polynomials modified for tetrahedra,
      // which are defined only on the reference/master tetrahedron.
      // Thus, to determine the high-order solution from the left and right
      // elements at the surface quadrature points, the basis functions
      // from the left and right elements are needed. For this, a
      // transformation to the reference coordinates is necessary, since
      // the basis functions are defined on the reference tetrahedron only.
      // Ref: [1] https://doi.org/10.1007/BF01060030
      //      [2] https://doi.org/10.1093/imamat/hxh111

      real detT_gp = 0.0;

      // transformation of the physical coordinates of the quadrature point
      // to reference space for the left element to be able to compute
      // basis functions on the left element.
      detT_gp = Jacobian( p1_l, gp, p3_l, p4_l );
      auto xi_l = detT_gp / detT_l;
      detT_gp = Jacobian( p1_l, p2_l, gp, p4_l );
      auto eta_l = detT_gp / detT_l;
      detT_gp = Jacobian( p1_l, p2_l, p3_l, gp );
      auto zeta_l = detT_gp / detT_l;

      // basis functions at igp for the left element
      auto B2l = 2.0 * xi_l + eta_l + zeta_l - 1.0;
      auto B3l = 3.0 * eta_l + zeta_l - 1.0;
      auto B4l = 4.0 * zeta_l - 1.0;

      // transformation of the physical coordinates of the quadrature point
      // to reference space for the right element
      detT_gp = Jacobian( p1_r, gp, p3_r, p4_r );
      auto xi_r = detT_gp / detT_r;
      detT_gp = Jacobian( p1_r, p2_r, gp, p4_r );
      auto eta_r = detT_gp / detT_r;
      detT_gp = Jacobian( p1_r, p2_r, p3_r, gp );
      auto zeta_r = detT_gp / detT_r;

      // basis functions at igp for the right element
      auto B2r = 2.0 * xi_r + eta_r + zeta_r - 1.0;
      auto B3r = 3.0 * eta_r + zeta_r - 1.0;
      auto B4r = 4.0 * zeta_r - 1.0;

      auto wt = wgp[igp] * geoFace(f,0,0);

      std::array< std::vector< real >, 2 >
        state{{ std::vector< real >( ncomp, 0.0 ),
                std::vector< real >( ncomp, 0.0 ) }};

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        auto lmark = c*(ndof-1);
        state[0][c] = U( el, mark, offset )
                      + limFunc(el, lmark+0, 0) * U(el, mark+1, offset) * B2l
                      + limFunc(el, lmark+1, 0) * U(el, mark+2, offset) * B3l
                      + limFunc(el, lmark+2, 0) * U(el, mark+3, offset) * B4l;
        state[1][c] = U( er, mark, offset )
                      + limFunc(er, lmark+0, 0) * U(er, mark+1, offset) * B2r
                      + limFunc(er, lmark+1, 0) * U(er, mark+2, offset) * B3r
                      + limFunc(er, lmark+2, 0) * U(er, mark+3, offset) * B4r;
      }

      // evaluate prescribed velocity (if any)
      auto v = vel( system, ncomp, gp[0], gp[1], gp[2] );
      // compute flux
      auto fl =
         flux( {{geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0)}}, state, v );

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;

        R(el, mark,   offset) -= wt * fl[c];
        R(el, mark+1, offset) -= wt * fl[c] * B2l;
        R(el, mark+2, offset) -= wt * fl[c] * B3l;
        R(el, mark+3, offset) -= wt * fl[c] * B4l;

        R(er, mark,   offset) += wt * fl[c];
        R(er, mark+1, offset) += wt * fl[c] * B2r;
        R(er, mark+2, offset) += wt * fl[c] * B3r;
        R(er, mark+3, offset) += wt * fl[c] * B4r;
      }
    }
  }
}
