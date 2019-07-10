// *****************************************************************************
/*!
  \file      src/PDE/Reconstruction.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Reconstruction for reconstructed discontinuous Galerkin methods
  \details   This file contains functions that reconstruct an "n"th order
    polynomial to an "n+1"th order polynomial using a least-squares
    reconstruction procedure.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "Vector.hpp"
#include "Reconstruction.hpp"

void
tk::intLeastSq_P0P1( ncomp_t ncomp,
                     ncomp_t offset,
                     const std::size_t rdof,
                     const inciter::FaceData& fd,
                     const Fields& geoElem,
                     const Fields& U,
                     std::vector< std::array< std::array< real, 3 >, 3 > >& lhs_ls,
                     std::vector< std::vector< std::array< real, 3 > > >& rhs_ls )
// *****************************************************************************
//  Compute internal surface contributions to the least-squares reconstruction
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoElem Element geometry array
//! \param[in] U Solution vector at recent time step
//! \param[in,out] lhs_ls LHS reconstruction matrix
//! \param[in,out] rhs_ls RHS reconstruction vector
// *****************************************************************************
{
  const auto& esuf = fd.Esuf();

  // Compute internal face contributions
  for (auto f=fd.Nbfac(); f<esuf.size()/2; ++f)
  {
    Assert( esuf[2*f] > -1 && esuf[2*f+1] > -1, "Interior element detected "
            "as -1" );

    std::size_t el = static_cast< std::size_t >(esuf[2*f]);
    std::size_t er = static_cast< std::size_t >(esuf[2*f+1]);

    // get a 3x3 system by applying the normal equation approach to the
    // least-squares overdetermined system
    std::array< real, 3 > wdeltax{{ 0.0, 0.0, 0.0 }};

    for (std::size_t idir=0; idir<3; ++idir)
      wdeltax[idir] = geoElem(er,idir+1,0)-geoElem(el,idir+1,0);

    for (std::size_t idir=0; idir<3; ++idir)
    {
      // rhs vector
      for (ncomp_t c=0; c<ncomp; ++c)
      {
        auto mark = c*rdof;
        rhs_ls[el][c][idir] +=
          wdeltax[idir] * (U(er,mark,offset)-U(el,mark,offset));
        rhs_ls[er][c][idir] +=
          wdeltax[idir] * (U(er,mark,offset)-U(el,mark,offset));
      }

      // lhs matrix
      for (std::size_t jdir=0; jdir<3; ++jdir)
      {
        lhs_ls[el][idir][jdir] += wdeltax[idir] * wdeltax[jdir];
        lhs_ls[er][idir][jdir] += wdeltax[idir] * wdeltax[jdir];
      }
    }
  }
}

void
tk::bndLeastSq_P0P1( ncomp_t system,
                     ncomp_t ncomp,
                     ncomp_t offset,
                     const std::size_t rdof,
                     const std::vector< bcconf_t >& bcconfig,
                     const inciter::FaceData& fd,
                     const Fields& geoFace,
                     const Fields& geoElem,
                     real t,
                     const StateFn& state,
                     const Fields& U,
                     std::vector< std::array< std::array< real, 3 >, 3 > >& lhs_ls,
                     std::vector< std::vector< std::array< real, 3 > > >& rhs_ls )
// *****************************************************************************
//  Compute boundary face contributions to the least-squares reconstruction
//! \param[in] system Equation system index
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] bcconfig BC configuration vector for multiple side sets
//! \param[in] fd Face connectivity and boundary conditions object
//! \param[in] geoFace Face geometry array
//! \param[in] geoElem Element geometry array
//! \param[in] t Physical time
//! \param[in] state Function to evaluate the left and right solution state at
//!   boundaries
//! \param[in] U Solution vector at recent time step
//! \param[in,out] lhs_ls LHS reconstruction matrix
//! \param[in,out] rhs_ls RHS reconstruction vector
// *****************************************************************************
{
  IGNORE(rdof);
  const auto& bface = fd.Bface();
  const auto& esuf = fd.Esuf();

  for (const auto& s : bcconfig) {       // for all bc sidesets
    auto bc = bface.find( std::stoi(s) );// faces for side set
    if (bc != end(bface))
    {
      // Compute boundary face contributions
      for (const auto& f : bc->second)
      {
        Assert( esuf[2*f+1] == -1, "outside boundary element not -1" );

        std::size_t el = static_cast< std::size_t >(esuf[2*f]);

        // arrays for quadrature points
        std::array< real, 3 >
          fc{{ geoFace(f,4,0), geoFace(f,5,0), geoFace(f,6,0)  }};
        std::array< real, 3 >
          fn{{ geoFace(f,1,0), geoFace(f,2,0), geoFace(f,3,0) }};

        // Compute the state variables at the left element
        std::vector< real >B(1,1.0);
        auto ul = eval_state( ncomp, offset, 1, 1, el, U, B );

        Assert( ul.size() == ncomp, "Size mismatch" );

        // Compute the state at the face-center using BC
        auto ustate = state( system, ncomp, ul, fc[0], fc[1], fc[2], t, fn );

        std::array< real, 3 > wdeltax{{ 0.0, 0.0, 0.0 }};

        for (std::size_t idir=0; idir<3; ++idir)
          wdeltax[idir] = fc[idir]-geoElem(el,1+idir,0);

        for (std::size_t idir=0; idir<3; ++idir)
        {
          // rhs vector
          for (ncomp_t c=0; c<ncomp; ++c) 
            rhs_ls[el][c][idir] +=
              wdeltax[idir] * (ustate[1][c]-ustate[0][c]);

          // lhs matrix
          for (std::size_t jdir=0; jdir<3; ++jdir)
            lhs_ls[el][idir][jdir] += wdeltax[idir] * wdeltax[jdir];
        }
      }
    }
  }
}

void
tk::solveLeastSq_P0P1( ncomp_t ncomp,
  ncomp_t offset,
  const std::size_t rdof,
  const std::vector< std::array< std::array< real, 3 >, 3 > >& lhs,
  const std::vector< std::vector< std::array< real, 3 > > >& rhs,
  Fields& U )
// *****************************************************************************
//  Solve 3x3 system for least-squares reconstruction
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Offset this PDE system operates from
//! \param[in] rdof Maximum number of reconstructed degrees of freedom
//! \param[in] lhs LHS reconstruction matrix
//! \param[in] rhs RHS reconstruction vector
//! \param[in,out] U Solution vector at recent time step
// *****************************************************************************
{
  for (std::size_t e=0; e<lhs.size(); ++e)
  {
    auto de = tk::determinant3by3( lhs[e] );

    for (ncomp_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      real nu(0.0);

      // solve system using Cramer's rule

      nu = tk::determinant3by3( {{{{rhs[e][c][0], lhs[e][0][1], lhs[e][0][2]}},
                                  {{rhs[e][c][1], lhs[e][1][1], lhs[e][1][2]}},
                                  {{rhs[e][c][2], lhs[e][2][1], lhs[e][2][2]}}}} );
      U(e,mark+1,offset) = nu/de;

      nu = tk::determinant3by3( {{{{lhs[e][0][0], rhs[e][c][0], lhs[e][0][2]}},
                                  {{lhs[e][1][0], rhs[e][c][1], lhs[e][1][2]}},
                                  {{lhs[e][2][0], rhs[e][c][2], lhs[e][2][2]}}}} );
      U(e,mark+2,offset) = nu/de;

      nu = tk::determinant3by3( {{{{lhs[e][0][0], lhs[e][0][1], rhs[e][c][0]}},
                                  {{lhs[e][1][0], lhs[e][1][1], rhs[e][c][1]}},
                                  {{lhs[e][2][0], lhs[e][2][1], rhs[e][c][2]}}}} );
      U(e,mark+3,offset) = nu/de;
    }
  }
}

void
tk::transform_P0P1( ncomp_t ncomp,
                    ncomp_t offset,
                    std::size_t rdof,
                    std::size_t nelem,
                    const std::vector< std::size_t >& inpoel,
                    const UnsMesh::Coords& coord,
                    Fields& U )
// *****************************************************************************
//  Transform the reconstructed P1-derivatives to the Dubiner dofs
//! \param[in] ncomp Number of scalar components in this PDE system
//! \param[in] offset Index for equation systems
//! \param[in] rdof Total number of reconstructed dofs
//! \param[in] nelem Total number of elements
//! \param[in] inpoel Element-node connectivity
//! \param[in] coord Array of nodal coordinates
//! \param[in,out] U Second-order solution vector which gets transformed to
//!   the Dubiner reference space
// *****************************************************************************
{
  const auto& cx = coord[0];
  const auto& cy = coord[1];
  const auto& cz = coord[2];

  for (std::size_t e=0; e<nelem; ++e)
  {
    // Extract the element coordinates
    std::array< std::array< real, 3>, 4 > coordel {{
      {{ cx[ inpoel[4*e  ] ], cy[ inpoel[4*e  ] ], cz[ inpoel[4*e  ] ] }},
      {{ cx[ inpoel[4*e+1] ], cy[ inpoel[4*e+1] ], cz[ inpoel[4*e+1] ] }},
      {{ cx[ inpoel[4*e+2] ], cy[ inpoel[4*e+2] ], cz[ inpoel[4*e+2] ] }},
      {{ cx[ inpoel[4*e+3] ], cy[ inpoel[4*e+3] ], cz[ inpoel[4*e+3] ] }}
    }};

    auto jacInv =
      tk::inverseJacobian( coordel[0], coordel[1], coordel[2], coordel[3] );

    // Compute the derivatives of basis function for DG(P1)
    auto dBdx = tk::eval_dBdx_p1( rdof, jacInv );

    std::array< std::array< real, 3 >, 3> dBdxa;

    for (std::size_t i=0; i<3; ++i)
      for (std::size_t j=0; j<3; ++j)
        dBdxa[i][j] = dBdx[i][j+1];

    auto de = tk::determinant3by3( dBdxa );

    for (ncomp_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      std::array< real, 3 > rhs{{ U(e,mark+1,offset),
                                  U(e,mark+2,offset),
                                  U(e,mark+3,offset) }};

      // solve system using Cramer's rule
      real nu(0.0);

      nu = tk::determinant3by3( {{{{rhs[0], dBdxa[0][1], dBdxa[0][2]}},
                                  {{rhs[1], dBdxa[1][1], dBdxa[1][2]}},
                                  {{rhs[2], dBdxa[2][1], dBdxa[2][2]}}}} );
      auto ux = nu/de;

      nu = tk::determinant3by3( {{{{dBdxa[0][0], rhs[0], dBdxa[0][2]}},
                                  {{dBdxa[1][0], rhs[1], dBdxa[1][2]}},
                                  {{dBdxa[2][0], rhs[2], dBdxa[2][2]}}}} );
      auto uy = nu/de;

      nu = tk::determinant3by3( {{{{dBdxa[0][0], dBdxa[0][1], rhs[0]}},
                                  {{dBdxa[1][0], dBdxa[1][1], rhs[1]}},
                                  {{dBdxa[2][0], dBdxa[2][1], rhs[2]}}}} );
      auto uz = nu/de;

      // replace physical derivatives with transformed dofs
      U(e,mark+1,offset) = ux;
      U(e,mark+2,offset) = uy;
      U(e,mark+3,offset) = uz;
    }
  }
}
