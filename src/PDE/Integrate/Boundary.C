// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Surface.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Functions for computing boundary surface integrals of a system
     of PDEs in DG methods
  \details   This file contains functionality for computing boundary surface
     integrals of a system of PDEs used in discontinuous Galerkin methods for
     various orders of numerical representation.
*/
// *****************************************************************************

#include <array>

#include "Boundary.h"
#include "Vector.h"
#include "Quadrature.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

void
tk::bndSurfIntP0( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const std::vector< std::size_t >& faces,
                  const std::vector< int >& esuf,
                  const Fields& geoFace,
                  real t,
                  const RiemannFluxFn& flux,
                  const VelFn& vel,
                  const StateFn& state,
                  const Fields& U,
                  Fields& R )
// *****************************************************************************
//  Compute boundary surface integral for a number of faces for DG(P0)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] faces Face IDs at which to compute surface integral
//! \param[in] esuf Elements surrounding face, see tk::genEsuf()
//! \param[in] geoFace Face geometry array
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  for (const auto& f : faces) {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );
    auto farea = geoFace(f,0,0);

    std::vector< real > ugp;
    for (ncomp_t c=0; c<ncomp; ++c)
      ugp.push_back( U(el, c, offset) );

    std::array< real, 3 >
      fn{{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }};
    auto xc = geoFace(f,4,0);
    auto yc = geoFace(f,5,0);
    auto zc = geoFace(f,6,0);

    auto fl = flux( fn,
                    state( system, ncomp, ugp, xc, yc, zc, t, fn ),
                    vel( xc, yc, zc, system, ncomp ) );

    for (ncomp_t c=0; c<ncomp; ++c)
      R(el, c, offset) -= farea * fl[c];
  }
}

void
tk::sidesetIntP0( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const std::vector< bcconf_t >& bcconfig,
                  const std::map< int, std::vector< std::size_t > >& bface,
                  const std::vector< int >& esuf,
                  const Fields& geoFace,
                  real t,
                  const RiemannFluxFn& flux,
                  const VelFn& vel,
                  const StateFn& state,
                  const Fields& U,
                  Fields& R )
// *****************************************************************************
//  Compute boundary surface flux integrals for a given boundary type for DG(P0)
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] bface Boundary faces side-set information
//! \param[in] esuf Elements surrounding faces
//! \param[in] geoFace Face geometry array
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in,out] R Right-hand side vector added to
// *****************************************************************************
{
  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find( std::stoi(s) );// faces for side set
    if (bc != end(bface))
      bndSurfIntP0( system, ncomp, offset, bc->second, esuf, geoFace, t, flux,
                    vel, state, U, R );
  }
}

// *****************************************************************************
//  Compute boundary surface integral for a number of faces for DG(P1)
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] faces Face IDs at which to compute surface integral
//! \param[in] esuf Elements surrounding face, see tk::genEsuf()
//! \param[in] geoFace Face geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] inpofa Face-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in,out] R Right-hand side vector computed
//! \tparam State Policy class providing the left and right state at
//!   boundaries by its member function State::LR()
// *****************************************************************************
void
tk::bndSurfIntP1( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const std::vector< std::size_t >& faces,
                  const std::vector< int >& esuf,
                  const Fields& geoFace,
                  const std::vector< std::size_t >& inpoel,
                  const std::vector< std::size_t >& inpofa,
                  const UnsMesh::Coords& coord,
                  real t,
                  const RiemannFluxFn& flux,
                  const VelFn& vel,
                  const StateFn& state,
                  const Fields& U,
                  const Fields& limFunc,
                  Fields& R )
{
  // Number of integration points
  constexpr std::size_t NG = 3;

  // arrays for quadrature points
  std::array< std::array< real, NG >, 2 > coordgp;
  std::array< real, NG > wgp;

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();

  // get quadrature point weights and coordinates for triangle
  GaussQuadratureTri( coordgp, wgp );

  for (const auto& f : faces) {
    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

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

    auto detT_l = Jacobian( p1_l, p2_l, p3_l, p4_l );

    auto x1 = cx[ inpofa[3*f]   ];
    auto y1 = cy[ inpofa[3*f]   ];
    auto z1 = cz[ inpofa[3*f]   ];

    auto x2 = cx[ inpofa[3*f+1] ];
    auto y2 = cy[ inpofa[3*f+1] ];
    auto z2 = cz[ inpofa[3*f+1] ];

    auto x3 = cx[ inpofa[3*f+2] ];
    auto y3 = cy[ inpofa[3*f+2] ];
    auto z3 = cz[ inpofa[3*f+2] ];

     std::array< real, 3 >
      fn{{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }};

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
      auto xgp = x1*shp1 + x2*shp2 + x3*shp3;
      auto ygp = y1*shp1 + y2*shp2 + y3*shp3;
      auto zgp = z1*shp1 + z2*shp2 + z3*shp3;

      real detT_gp(0);

      // transformation of the physical coordinates of the quadrature point
      // to reference space for the left element to be able to compute
      // basis functions on the left element.
      detT_gp = Jacobian( p1_l, {{ xgp, ygp, zgp }}, p3_l, p4_l );
      auto xi_l = detT_gp / detT_l;
      detT_gp = Jacobian( p1_l, p2_l, {{ xgp, ygp, zgp }}, p4_l );
      auto eta_l = detT_gp / detT_l;
      detT_gp = Jacobian( p1_l, p2_l, p3_l, {{ xgp, ygp, zgp }} );
      auto zeta_l = detT_gp / detT_l;

      // basis functions at igp for the left element
      auto B2l = 2.0 * xi_l + eta_l + zeta_l - 1.0;
      auto B3l = 3.0 * eta_l + zeta_l - 1.0;
      auto B4l = 4.0 * zeta_l - 1.0;

      auto wt = wgp[igp] * geoFace(f,0,0);

      std::vector< real > ugp( ncomp );

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        auto lmark = c*(ndof-1);
        ugp[c] = U(el, mark,   offset)
                 + limFunc(el, lmark+0, 0) * U(el, mark+1, offset) * B2l
                 + limFunc(el, lmark+1, 0) * U(el, mark+2, offset) * B3l
                 + limFunc(el, lmark+2, 0) * U(el, mark+3, offset) * B4l;
      }

      auto fl = flux( fn,
                      state( system, ncomp, ugp, xgp, ygp, zgp, t, fn ),
                      vel( xgp, ygp, zgp, system, ncomp ) );

      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*ndof;
        R(el, mark,   offset) -= wt * fl[c];
        R(el, mark+1, offset) -= wt * fl[c] * B2l;
        R(el, mark+2, offset) -= wt * fl[c] * B3l;
        R(el, mark+3, offset) -= wt * fl[c] * B4l;
      }
    }
  }
}

void
tk::sidesetIntP1( ncomp_t system,
                  ncomp_t ncomp,
                  ncomp_t offset,
                  const std::vector< bcconf_t >& bcconfig,
                  const std::map< int, std::vector< std::size_t > >& bface,
                  const std::vector< int >& esuf,
                  const Fields& geoFace,
                  const std::vector< std::size_t >& inpoel,
                  const std::vector< std::size_t >& inpofa,
                  const UnsMesh::Coords& coord,
                  real t,
                  const RiemannFluxFn& flux,
                  const VelFn& vel,
                  const StateFn& state,
                  const Fields& U,
                  const Fields& limFunc,
                  Fields& R )
// *****************************************************************************
//! Compute boundary surface flux integrals for a given boundary type for DG(P1)
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] bface Boundary faces side-set information
//! \param[in] esuf Elements surrounding faces
//! \param[in] geoFace Face geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] inpofa Face-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] limFunc Limiter function for higher-order solution dofs
//! \param[in,out] R Right-hand side vector computed
// *****************************************************************************
{
  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find( std::stoi(s) );// faces for side set
    if (bc != end(bface))
      bndSurfIntP1( system, ncomp, offset, bc->second, esuf, geoFace, inpoel,
                    inpofa, coord, t, flux, vel, state, U, limFunc, R );
  }
}
