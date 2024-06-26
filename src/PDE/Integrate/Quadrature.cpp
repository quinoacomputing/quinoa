// *****************************************************************************
/*!
  \file      src/PDE/Integrate/Quadrature.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Quadrature coordinates and weights for numerical integration
  \details   This file contains functions that provide Gauss quadrature
    coordinates and weights for numerical integration.
*/
// *****************************************************************************

#include "Quadrature.hpp"

void
tk::GaussQuadratureTet( const std::size_t NG,
                        std::array< std::vector< real >, 3>& coordgp,
                        std::vector< real >& wgp )
// *****************************************************************************
//! Initialize Gaussian quadrature points locations and weights for a tetrahedron
//! \param[in] NG number of quadrature points
//! \param[in,out] coordgp 3 spatial coordinates of quadrature points
//! \param[in,out] wgp Weights of quadrature points
// *****************************************************************************
{
  Assert( coordgp[0].size() == NG, "Size mismatch" );
  Assert( coordgp[1].size() == NG, "Size mismatch" );
  Assert( coordgp[1].size() == NG, "Size mismatch" );
  Assert( wgp.size() == NG, "Size mismatch" );

  switch( NG )
  {
    case 1:
    {
      coordgp[0][0] = 0.25;
      coordgp[1][0] = 0.25;
      coordgp[2][0] = 0.25;
      wgp[0]        = 1.0;
    }
    break;

    case 4:
    {
      const real a1 = 0.5854101966249685;
      const real a2 = 0.1381966011250105;

      coordgp[0][0] = a2;
      coordgp[1][0] = a2;
      coordgp[2][0] = a2;
      wgp[0]        = 0.25;

      coordgp[0][1] = a1;
      coordgp[1][1] = a2;
      coordgp[2][1] = a2;
      wgp[1]        = 0.25;

      coordgp[0][2] = a2;
      coordgp[1][2] = a1;
      coordgp[2][2] = a2;
      wgp[2]        = 0.25;

      coordgp[0][3] = a2;
      coordgp[1][3] = a2;
      coordgp[2][3] = a1;
      wgp[3]        = 0.25;
    }
    break;

    case 5:
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
    break;

    case 11:
    {
      const tk::real c1 = 0.3994035761667992;
      const tk::real c2 = 0.1005964238332008;
      const tk::real c3 = 343.0 / 7500.0;
      const tk::real c4 = 56.0 / 375.0;

      coordgp[0][0] = 0.25;
      coordgp[1][0] = 0.25;
      coordgp[2][0] = 0.25;
      wgp[0]        = -148.0/1875.0;

      coordgp[0][1] = 11.0/14.0;
      coordgp[1][1] = 1.0/14.0;
      coordgp[2][1] = 1.0/14.0;
      wgp[1]        = c3;

      coordgp[0][2] = 1.0/14.0;
      coordgp[1][2] = 11.0/14.0;
      coordgp[2][2] = 1.0/14.0;
      wgp[2]        = c3;

      coordgp[0][3] = 1.0/14.0;
      coordgp[1][3] = 1.0/14.0;
      coordgp[2][3] = 11.0/14.0;
      wgp[3]        = c3;

      coordgp[0][4] = 1.0/14.0;
      coordgp[1][4] = 1.0/14.0;
      coordgp[2][4] = 1.0/14.0;
      wgp[4]        = c3;

      coordgp[0][5] = c1;
      coordgp[1][5] = c1;
      coordgp[2][5] = c2;
      wgp[5]        = c4;

      coordgp[0][6] = c1;
      coordgp[1][6] = c2;
      coordgp[2][6] = c1;
      wgp[6]        = c4;

      coordgp[0][7] = c1;
      coordgp[1][7] = c2;
      coordgp[2][7] = c2;
      wgp[7]        = c4;

      coordgp[0][8] = c2;
      coordgp[1][8] = c1;
      coordgp[2][8] = c1;
      wgp[8]        = c4;

      coordgp[0][9] = c2;
      coordgp[1][9] = c1;
      coordgp[2][9] = c2;
      wgp[9]        = c4;

      coordgp[0][10] = c2;
      coordgp[1][10] = c2;
      coordgp[2][10] = c1;
      wgp[10]        = c4;
    }
    break;

    case 14:
    {
      const tk::real a = 0.0673422422100983;
      const tk::real b = 0.3108859192633005;
      const tk::real c = 0.7217942490673264;
      const tk::real d = 0.0927352503108912;
      const tk::real e = 0.4544962958743506;
      const tk::real f = 0.0455037041256494;
      const tk::real p = 0.1126879257180162;
      const tk::real q = 0.0734930431163619;
      const tk::real r = 0.0425460207770812;

      coordgp[0][0] = a;
      coordgp[1][0] = b;
      coordgp[2][0] = b;
      wgp[0]        = p;

      coordgp[0][1] = b;
      coordgp[1][1] = a;
      coordgp[2][1] = b;
      wgp[1]        = p;

      coordgp[0][2] = b;
      coordgp[1][2] = b;
      coordgp[2][2] = a;
      wgp[2]        = p;

      coordgp[0][3] = b;
      coordgp[1][3] = b;
      coordgp[2][3] = b;
      wgp[3]        = p;

      coordgp[0][4] = c;
      coordgp[1][4] = d;
      coordgp[2][4] = d;
      wgp[4]        = q;

      coordgp[0][5] = d;
      coordgp[1][5] = c;
      coordgp[2][5] = d;
      wgp[5]        = q;

      coordgp[0][6] = d;
      coordgp[1][6] = d;
      coordgp[2][6] = c;
      wgp[6]        = q;

      coordgp[0][7] = d;
      coordgp[1][7] = d;
      coordgp[2][7] = d;
      wgp[7]        = q;

      coordgp[0][8] = e;
      coordgp[1][8] = e;
      coordgp[2][8] = f;
      wgp[8]        = r;

      coordgp[0][9] = e;
      coordgp[1][9] = f;
      coordgp[2][9] = e;
      wgp[9]        = r;

      coordgp[0][10] = e;
      coordgp[1][10] = f;
      coordgp[2][10] = f;
      wgp[10]        = r;

      coordgp[0][11] = f;
      coordgp[1][11] = e;
      coordgp[2][11] = e;
      wgp[11]        = r;

      coordgp[0][12] = f;
      coordgp[1][12] = e;
      coordgp[2][12] = f;
      wgp[12]        = r;

      coordgp[0][13] = f;
      coordgp[1][13] = f;
      coordgp[2][13] = e;
      wgp[13]        = r;
    }
    break;
  }
}

void
tk::GaussQuadratureTri( const std::size_t NG,
                        std::array< std::vector< real >, 2>& coordgp,
                        std::vector< real >& wgp )
// *****************************************************************************
//! Initialize Gaussian quadrature points locations and weights for a triangle
//! \param[in] NG number of quadrature points
//! \param[in,out] coordgp 2 spatial coordinates of quadrature points
//! \param[in,out] wgp Weights of quadrature points
// *****************************************************************************
{
  Assert( coordgp[0].size() == NG, "Size mismatch" );
  Assert( coordgp[1].size() == NG, "Size mismatch" );
  Assert( wgp.size() == NG, "Size mismatch" );
  
  switch( NG )
  {
    case 1:
      coordgp[0][0] = 1.0/3.0;
      coordgp[1][0] = 1.0/3.0;
      wgp[0]        = 1.0;
      break;

    case 3:
      coordgp[0][0] = 2.0/3.0;
      coordgp[1][0] = 1.0/6.0;
      wgp[0]        = 1.0/3.0;

      coordgp[0][1] = 1.0/6.0;
      coordgp[1][1] = 2.0/3.0;
      wgp[1]        = 1.0/3.0;

      coordgp[0][2] = 1.0/6.0;
      coordgp[1][2] = 1.0/6.0;
      wgp[2]        = 1.0/3.0;
      break;

    case 4:
      coordgp[0][0] = 1.0/3.0;
      coordgp[1][0] = 1.0/3.0;
      wgp[0]        = -27.0/48.0;

      coordgp[0][1] = 1.0/5.0;
      coordgp[1][1] = 1.0/5.0;
      wgp[1]        = 25.0/48.0;

      coordgp[0][2] = 3.0/5.0;
      coordgp[1][2] = 1.0/5.0;
      wgp[2]        = 25.0/48.0;

      coordgp[0][3] = 1.0/5.0;
      coordgp[1][3] = 3.0/5.0;
      wgp[3]        = 25.0/48.0;
      break;

    case 6:
      const tk::real c1 = 0.816847572980459;
      const tk::real c2 = 0.091576213509771;
      const tk::real c3 = 0.091576213509771;
      const tk::real c4 = 0.108103018168070;
      const tk::real c5 = 0.445948490915965;
      const tk::real c6 = 0.445948490915965;
      const tk::real w1 = 0.109951743655322;
      const tk::real w2 = 0.223381589678011;

      coordgp[0][0] = c1;
      coordgp[1][0] = c2;
      wgp[0]        = w1;

      coordgp[0][1] = c2;
      coordgp[1][1] = c3;
      wgp[1]        = w1;

      coordgp[0][2] = c3;
      coordgp[1][2] = c1;
      wgp[2]        = w1;

      coordgp[0][3] = c4;
      coordgp[1][3] = c5;
      wgp[3]        = w2;

      coordgp[0][4] = c5;
      coordgp[1][4] = c6;
      wgp[4]        = w2;

      coordgp[0][5] = c6;
      coordgp[1][5] = c4;
      wgp[5]        = w2;
      break;
  }
}
