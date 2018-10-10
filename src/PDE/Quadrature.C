// *****************************************************************************
/*!
  \file      src/PDE/Quadrature.C
  \copyright 2016-2018, Los Alamos National Security, LLC.
  \brief     Quadrature coordinates and weights for numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for numerical integration.
*/
// *****************************************************************************

#include <array>
#include <vector>

#include "Quadrature.h"

void
inciter::GaussQuadratureTet(
  std::array< std::array< tk::real, 5 >, 3 >& coordgp,
  std::array< tk::real, 5 >& wgp )
// *****************************************************************************
//  Five Gaussian quadrature points locations and weights for a tetrahedron
//! \param[in,out] coordgp 3 spatial coordinates of 5 quadrature points
//! \param[in,out] wgp 5 weights of quadrature points
// *****************************************************************************
{
  coordgp[0][0] = 0.25;
  coordgp[1][0] = 0.25;
  coordgp[2][0] = 0.25;
  wgp[0]        = -12.0/15.0;

  coordgp[0][1] = 1.0/6.0;
  coordgp[1][1] = 1.0/6.0;
  coordgp[2][1] = 1.0/6.0;
  wgp[1]        = 9.0/20.0;

  coordgp[0][2] = 0.5;
  coordgp[1][2] = 1.0/6.0;
  coordgp[2][2] = 1.0/6.0;
  wgp[2]        = 9.0/20.0;

  coordgp[0][3] = 1.0/6.0;
  coordgp[1][3] = 0.5;
  coordgp[2][3] = 1.0/6.0;
  wgp[3]        = 9.0/20.0;

  coordgp[0][4] = 1.0/6.0;
  coordgp[1][4] = 1.0/6.0;
  coordgp[2][4] = 0.5;
  wgp[4]        = 9.0/20.0;
}

void
inciter::GaussQuadratureTet(
  std::array< std::array< tk::real, 4 >, 3 >& coordgp,
  std::array< tk::real, 4 >& wgp )
// *****************************************************************************
//  Five Gaussian quadrature points locations and weights for a tetrahedron
//! \param[in,out] coordgp 3 spatial coordinates of 5 quadrature points
//! \param[in,out] wgp 5 weights of quadrature points
// *****************************************************************************
{
  const tk::real c1 = 0.5854101966249685;
  const tk::real c2 = 0.1381966011250105;

  coordgp[0][0] = c2;
  coordgp[1][0] = c2;
  coordgp[2][0] = c2;
  wgp[0]        = 0.25;

  coordgp[0][1] = c1;
  coordgp[1][1] = c2;
  coordgp[2][1] = c2;
  wgp[1]        = 0.25;

  coordgp[0][2] = c2;
  coordgp[1][2] = c1;
  coordgp[2][2] = c2;
  wgp[2]        = 0.25;

  coordgp[0][3] = c2;
  coordgp[1][3] = c2;
  coordgp[2][3] = c1;
  wgp[3]        = 0.25;
}

void
inciter::GaussQuadratureTri(
  std::array< std::array< tk::real, 3 >, 2 >& coordgp,
  std::array< tk::real, 3 >& wgp )
// *****************************************************************************
//! Three Gaussian quadrature points locations and weights for a triangle
//! \param[in,out] coordgp 2 spatial coordinates of 3 quadrature points
//! \param[in,out] wgp 3 weights of quadrature points
// *****************************************************************************
{
  coordgp[0][0] = 2.0/3.0;
  coordgp[1][0] = 1.0/6.0;
  wgp[0]        = 1.0/3.0;

  coordgp[0][1] = 1.0/6.0;
  coordgp[1][1] = 2.0/3.0;
  wgp[1]        = 1.0/3.0;

  coordgp[0][2] = 1.0/6.0;
  coordgp[1][2] = 1.0/6.0;
  wgp[2]        = 1.0/3.0;
}
