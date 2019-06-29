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
#include "Inciter/InputDeck/InputDeck.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void
leastSquares_P0P1( const std::vector< int >& esuel,
                   inciter::ncomp_t offset,
                   const tk::Fields& geoElem,
                   tk::Fields& U )
// *****************************************************************************
//  Least-squares reconstruction for rDG(P0P1)
//! \param[in] esuel Elements surrounding elements
//! \param[in] offset Index for equation systems
//! \param[in] geoElem Element geometry array
//! \param[in,out] U First-order solution vector which gets reconstructed to
//!   second-order by the least-squares procedure
// *****************************************************************************
{
  const auto rdof = inciter::g_inputdeck.get< tag::discr, tag::rdof >();
  std::size_t ncomp = U.nprop()/rdof;

  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    for (inciter::ncomp_t c=0; c<ncomp; ++c)
    {
      auto mark = c*rdof;
      auto uc = U(e,mark,offset);
      std::array< tk::real, 3 > xe = {{ geoElem(e,1,0),
                                        geoElem(e,2,0),
                                        geoElem(e,3,0) }};


      // get a 3x3 system by applying the normal equation approach to the
      // least-squares overdetermined system
      std::array< tk::real, 3 > wdeltax{{ 0.0, 0.0, 0.0 }};
      std::array< tk::real, 3 > rhs_ls{{ 0.0, 0.0, 0.0 }};
      std::array< std::array< tk::real, 3 >, 3 > lhs_ls;
      for (std::size_t idir=0; idir<3; ++idir)
        lhs_ls[idir].fill(0.0);

      for (std::size_t j=0; j<4; ++j)
      {
        auto nel = esuel[4*e+j];
        if (nel == -1) continue;
        auto neL = static_cast< std::size_t >(nel);

        for (std::size_t idir=0; idir<3; ++idir)
          wdeltax[idir] = geoElem(neL,idir+1,0)-xe[idir];

        for (std::size_t idir=0; idir<3; ++idir)
        {
          rhs_ls[idir] += wdeltax[idir] * (U(neL,mark,offset)-uc);
          for (std::size_t jdir=0; jdir<3; ++jdir)
            lhs_ls[idir][jdir] += wdeltax[idir] * wdeltax[jdir];
        }
      }
      // solve system using Cramer's rule
      tk::real nu(0.0);
      auto de = tk::determinant3by3( lhs_ls );

      nu = tk::determinant3by3( {{{{rhs_ls[0], lhs_ls[0][1], lhs_ls[0][2]}},
                                  {{rhs_ls[1], lhs_ls[1][1], lhs_ls[1][2]}},
                                  {{rhs_ls[2], lhs_ls[2][1], lhs_ls[2][2]}}}} );
      U(e,mark+1,offset) = nu/de;

      nu = tk::determinant3by3( {{{{lhs_ls[0][0], rhs_ls[0], lhs_ls[0][2]}},
                                  {{lhs_ls[1][0], rhs_ls[1], lhs_ls[1][2]}},
                                  {{lhs_ls[2][0], rhs_ls[2], lhs_ls[2][2]}}}} );
      U(e,mark+2,offset) = nu/de;

      nu = tk::determinant3by3( {{{{lhs_ls[0][0], lhs_ls[0][1], rhs_ls[0]}},
                                  {{lhs_ls[1][0], lhs_ls[1][1], rhs_ls[1]}},
                                  {{lhs_ls[2][0], lhs_ls[2][1], rhs_ls[2]}}}} );
      U(e,mark+3,offset) = nu/de;
    }
  }
}

} // inciter::
