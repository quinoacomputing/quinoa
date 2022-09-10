// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Source.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Functions for computing integrals of an arbitrary source term of a
     for single-material compressible flow, CompFlow using DG methods
  \details   This file contains functionality for computing integrals of an
     arbitrary source term for single-material compressible flow, CompFlow with
     discontinuous Galerkin methods for various orders of numerical
     representation.
*/
// *****************************************************************************

#include <vector>

#include "Source.hpp"
#include "Quadrature.hpp"

void
tk::srcInt( ncomp_t system,
            const std::vector< inciter::EoS_Base* >& mat_blk,
            real t,
            const std::size_t ndof,
            const std::size_t nelem,
            const std::vector< std::size_t >& inpoel,
            const UnsMesh::Coords& coord,
            const Fields& geoElem,
            const SrcFn& src,
            const std::vector< std::size_t >& ndofel,
            Fields& R,
            std::size_t nmat )
// *****************************************************************************
//  Compute source term integrals for DG
//! \param[in] system Equation system index
//! \param[in] mat_blk Material block EOS
//! \param[in] t Physical time
//! \param[in] ndof Maximum number of degrees of freedom
//! \param[in] nelem Maximum number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in] geoElem Element geometry array
//! \param[in] src Source function to use
//! \param[in] ndofel Vector of local number of degrees of freedome
//! \param[in,out] R Right-hand side vector computed
//! \param[in] nmat Number of materials. A default is set to 1, so that calling
//!   code for single material systems primitive quantities does not need to
//!   specify this argument.
// *****************************************************************************
{
  auto ncomp = R.nprop()/ndof;
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

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
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }} }};

    for (std::size_t igp=0; igp<ng; ++igp)
    {
      // Compute the coordinates of quadrature point at physical domain
      auto gp = eval_gp( igp, coordel, coordgp );

      // Compute the basis function
      auto B =
        eval_basis( ndofel[e], coordgp[0][igp], coordgp[1][igp], coordgp[2][igp] );

      // Compute the source term variable
      std::vector< real > s(ncomp, 0.0);
      src( system, nmat, mat_blk, gp[0], gp[1], gp[2], t, s );

      auto wt = wgp[igp] * geoElem(e, 0);

      update_rhs( ndof, ndofel[e], wt, e, B, s, R );
    }
  }
}

void
tk::update_rhs( const std::size_t ndof,
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
  Assert( B.size() == ndof_el, "Size mismatch for basis function" );

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
      }
    }
  }
}

void
tk::srcIntFV( ncomp_t system,
              const std::vector< inciter::EoS_Base* >& mat_blk,
              real t,
              const std::size_t nelem,
              const Fields& geoElem,
              const SrcFn& src,
              Fields& R,
              std::size_t nmat )
// *****************************************************************************
//  Compute source term integrals for DG
//! \param[in] system Equation system index
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
    src( system, nmat, mat_blk, geoElem(e,1), geoElem(e,2), geoElem(e,3), t, s );

    // Add the source term to the rhs
    for (ncomp_t c=0; c<ncomp; ++c)
    {
      R(e, c) += geoElem(e,0) * s[c];
    }
  }
}
