// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Boundary.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
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
#include "Vector.hpp"
#include "Quadrature.hpp"

void
tk::bndSurfInt( ncomp_t system,
                ncomp_t ncomp,
                std::size_t nmat,
                ncomp_t offset,
                const std::size_t ndof,
                const std::vector< bcconf_t >& bcconfig,
                const inciter::FaceData& fd,
                const Fields& geoFace,
                const std::vector< std::size_t >& inpoel,
                const UnsMesh::Coords& coord,
                real t,
                const RiemannFluxFn& flux,
                const VelFn& vel,
                const StateFn& state,
                const Fields& U,
                const std::vector< std::size_t >& ndofel,
                Fields& R,
                std::vector< std::vector< tk::real > >& riemannDeriv )
// *****************************************************************************
//! Compute boundary surface flux integrals for a given boundary type for DG
//! \details This function computes contributions from surface integrals along
//!   all faces for a particular boundary condition type, configured by the state
//!   function
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] nmat Number of materials in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] t Physical time
//! \param[in] flux Riemann flux function to use
//! \param[in] vel Function to use to query prescribed velocity (if any)
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in] ndofel Vector of local number of degrees of freedom
//! \param[in,out] R Right-hand side vector computed
//! \param[in,out] riemannDeriv Derivatives of partial-pressures and velocities
//!   computed from the Riemann solver for use in the non-conservative terms.
//!   These derivatives are used only for multi-material hydro and unused for
//!   single-material compflow and linear transport.
// *****************************************************************************
{
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();
  const auto& inpofa = fd.Inpofa();

  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  Assert( (nmat==1 ? riemannDeriv.empty() : true), "Non-empty Riemann "
          "derivative vector for single material compflow" );

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find( std::stoi(s) );// faces for side set
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
          fn{{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }};

        // Gaussian quadrature
        for (std::size_t igp=0; igp<ng; ++igp)
        {
          // Compute the coordinates of quadrature point at physical domain
          auto gp = eval_gp( igp, coordfa, coordgp );

          //Compute the basis functions for the left element
          auto B_l = eval_basis( ndofel[el],
            Jacobian( coordel_l[0], gp, coordel_l[2], coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], gp, coordel_l[3] ) / detT_l,
            Jacobian( coordel_l[0], coordel_l[1], coordel_l[2], gp ) / detT_l );

          auto wt = wgp[igp] * geoFace(f,0,0);

          // Compute the state variables at the left element
          auto ugp = eval_state( ncomp, offset, ndof, ndofel[el], el, U, B_l );

          Assert( ugp.size() == ncomp, "Size mismatch" );

          // Compute the numerical flux
          auto fl = flux( fn,
                      state( system, ncomp, ugp, gp[0], gp[1], gp[2], t, fn ),
                      vel( system, ncomp, gp[0], gp[1], gp[2] ) );

          // Add the surface integration term to the rhs
          update_rhs_bc( ncomp, nmat, offset, ndof, ndofel[el], wt, fn, el, fl,
                         B_l, R, riemannDeriv );
        }
      }
    }
  }
}

void
tk::update_rhs_bc ( ncomp_t ncomp,
                    std::size_t nmat,
                    ncomp_t offset,
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
//! \param[in] offset Offset this PDE system operates from
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
  Assert( B_l.size() == ndof_l, "Size mismatch" );

  for (ncomp_t c=0; c<ncomp; ++c)
  {
    auto mark = c*ndof;
    R(el, mark, offset) -= wt * fl[c];

    if(ndof_l > 1)          //DG(P1)
    {
      R(el, mark+1, offset) -= wt * fl[c] * B_l[1];
      R(el, mark+2, offset) -= wt * fl[c] * B_l[2];
      R(el, mark+3, offset) -= wt * fl[c] * B_l[3];
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

    // Divergence of velocity
    riemannDeriv[3*nmat][el] += wt * fl[ncomp+nmat];
  }
}
