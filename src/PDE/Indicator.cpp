// *****************************************************************************
/*!
  \file      src/PDE/Indicator.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Adaptive indicators for p-adaptive discontiunous Galerkin methods
  \details   This file contains functions that provide adaptive indicator
    function calculations for marking the number of degree of freedom of each 
    element.
*/
// *****************************************************************************

#include "Indicator.hpp"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

void eval_ndof( const std::size_t nunk,
                const std::vector< int >& esuel,
                const tk::Fields& unk,
                std::vector< std::size_t >& ndofel )
// *****************************************************************************
//! Evaluate the adaptive indicator and mark the ndof for each element
//! \param[in] nunk Number of unknowns
//! \param[in] esuel Elements surrounding elements
//! \param[in] unk Array of unknowns
//! \param[in,out] ndofel Vector of local number of degrees of freedome
// *****************************************************************************
{
  const auto ndof = inciter::g_inputdeck.get< tag::discr, tag::ndof >();
  const auto ncomp = unk.nprop() / ndof;
  const auto ndofmax = inciter::g_inputdeck.get< tag::pref, tag::ndofmax >();

  // The array storing the adaptive indicator for each elements
  std::vector< tk::real > Ind(nunk, 0);

  tk::real ErrMaxp2(0), ErrMaxp1(0), ErrMinp2(1), ErrMinp1(1);

  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    if(ndofel[e] > 1)
    {
      auto ng = tk::NGvol(ndofel[e]);

      // arrays for quadrature points
      std::array< std::vector< tk::real >, 3 > coordgp;
      std::vector< tk::real > wgp;

      coordgp[0].resize( ng );
      coordgp[1].resize( ng );
      coordgp[2].resize( ng );
      wgp.resize( ng );

      tk::GaussQuadratureTet( ng, coordgp, wgp );

      // The adaptive indicator used here is the spectral decay indicator. This
      // indicator calculates the difference between the prohjections of the 
      // numerical solutions on finite element space of order p and p-1. For 
      // further references see:
      // https://doi.org/10.1016/j.jcp.2018.09.045

      tk::real dU(0), U(0);

      // Gaussian quadrature
      for (std::size_t igp=0; igp<ng; ++igp)
      {
        // Compute the basis function
        auto B = tk::eval_basis( ndofel[e], coordgp[0][igp], coordgp[1][igp],
                                 coordgp[2][igp] );

        auto state = tk::eval_state( ncomp, 0, ndof, ndofel[e], e, unk, B );

        U += wgp[igp] * state[0] * state[0];

        if(ndofel[e] > 4)
        {
           auto dU_p2 = unk(e, 4, 0) * B[4]
                      + unk(e, 5, 0) * B[5]
                      + unk(e, 6, 0) * B[6]
                      + unk(e, 7, 0) * B[7]
                      + unk(e, 8, 0) * B[8]
                      + unk(e, 9, 0) * B[9];
 
           dU += wgp[igp] * dU_p2 * dU_p2;
        }
        else
        {
           auto dU_p1 = unk(e, 1, 0) * B[1]
                      + unk(e, 2, 0) * B[2]
                      + unk(e, 3, 0) * B[3];
 
           dU += wgp[igp] * dU_p1 * dU_p1;
        }
      }

      //std::cout << "dU = " << dU << "\t" << "U = " << U << std::endl;
      Ind[e] = log10( dU / U );

      if(ndofel[e] > 4)
      {
        if(Ind[e] > ErrMaxp2)
          ErrMaxp2 = Ind[e];
        if(Ind[e] < ErrMinp2)
          ErrMinp2 = Ind[e];
      }
      else if(ndofel[e] == 4)
      {
        if(Ind[e] > ErrMaxp1)
          ErrMaxp1 = Ind[e];
        if(Ind[e] < ErrMinp1)
          ErrMinp1 = Ind[e];
      }
    }
  }

  //std::cout << "ErrMaxp1 = " << ErrMaxp1 << std::endl;
  //std::cout << "ErrMinp1 = " << ErrMinp1 << std::endl; 
  //std::cout << "ErrMaxp2 = " << ErrMaxp2 << std::endl;
  //std::cout << "ErrMinp2 = " << ErrMinp2 << std::endl;

  // Marke the ndof according to the adaptive indicator
  for (std::size_t e=0; e<esuel.size()/4; ++e)
  {
    if(ndofel[e] == 4)
    {
      if(Ind[e] > -6 && ndofmax > 4)      // Refinement
       ndofel[e] = 10;
      if(Ind[e] < -6.5)                     // Derefinement
        ndofel[e] = 1;
    }
    else if(ndofel[e] > 4)
    {
      if(Ind[e] < -8)                     // Derefinement
        ndofel[e] = 4;
    }
  }
}

} // inciter::
