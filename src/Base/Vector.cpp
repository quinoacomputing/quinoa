// *****************************************************************************
/*!
  \file      src/Base/Vector.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************

#include <array>
#include <cmath>

#include "Vector.hpp"
#include "Exception.hpp"

std::array< tk::real, 3 >
tk::cross( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 )
// *****************************************************************************
//  Compute the cross-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Cross-product
// *****************************************************************************
{
  return {{ v1[1]*v2[2] - v2[1]*v1[2],
            v1[2]*v2[0] - v2[2]*v1[0],
            v1[0]*v2[1] - v2[0]*v1[1] }};
}

std::array< tk::real, 3 >
tk::crossdiv( const std::array< real, 3 >& v1,
              const std::array< real, 3 >& v2,
              real j )
// *****************************************************************************
//  Compute the cross-product of two vectors divided by a scalar
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] j Scalar to divide each component by
//! \return Cross-product divided by scalar
// *****************************************************************************
{
  return {{ (v1[1]*v2[2] - v2[1]*v1[2]) / j,
            (v1[2]*v2[0] - v2[2]*v1[0]) / j,
            (v1[0]*v2[1] - v2[0]*v1[1]) / j }};
}

tk::real
tk::dot( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 )
// *****************************************************************************
//  Compute the dot-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Dot-product
// *****************************************************************************
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

tk::real
tk::length( const std::array< real, 3 >& v )
// *****************************************************************************
//  Compute length of a vector
//! \param[in] v vector
//! \return length
// *****************************************************************************
{
  return std::sqrt( dot(v,v) );
}

void
tk::unit( std::array< real, 3 >& v )
// *****************************************************************************
// Scale vector to unit length
//! \param[in,out] v Vector to normalize
// *****************************************************************************
{
  auto len = length( v );
  Assert( len > std::numeric_limits< tk::real >::epsilon(), "div by zero" );
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

tk::real
tk::triple( const std::array< real, 3 >& v1,
            const std::array< real, 3 >& v2,
            const std::array< real, 3 >& v3 )
// *****************************************************************************
//  Compute the triple-product of three vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] v3 3rd vector
//! \return Triple-product
// *****************************************************************************
{
  return dot( v1, cross(v2,v3) );
}

std::array< tk::real, 3 >
tk::rotateX( const std::array< real, 3 >& v, real angle )
// *****************************************************************************
//! Rotate vector about X axis by -45 degress
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
// *****************************************************************************
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ 1.0,         0.0,          0.0 }},
        {{ 0.0,   cos(angle), -sin(angle) }},
        {{ 0.0,   sin(angle),  cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

std::array< tk::real, 3 >
tk::rotateY( const std::array< real, 3 >& v, real angle )
// *****************************************************************************
//! Rotate vector about Y axis by -45 degress
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
// *****************************************************************************
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle),  0.0, sin(angle) }},
        {{ 0.0,         1.0,        0.0 }},
        {{ -sin(angle), 0.0, cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

std::array< tk::real, 3 >
tk::rotateZ( const std::array< real, 3 >& v, real angle )
// *****************************************************************************
//! Rotate vector about Z axis by -45 degress
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
// *****************************************************************************
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle), -sin(angle), 0.0 }},
        {{ sin(angle),  cos(angle), 0.0 }},
        {{ 0.0,         0.0,        1.0 }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

tk::real
tk::Jacobian( const std::array< real, 3 >& v1,
              const std::array< real, 3 >& v2,
              const std::array< real, 3 >& v3,
              const std::array< real, 3 >& v4 )
// *****************************************************************************
//  Compute the determinant of the Jacobian of a coordinate transformation over
//  a tetrahedron
//! \param[in] v1 (x,y,z) coordinates of 1st vertex of the tetrahedron
//! \param[in] v2 (x,y,z) coordinates of 2nd vertex of the tetrahedron
//! \param[in] v3 (x,y,z) coordinates of 3rd vertex of the tetrahedron
//! \param[in] v4 (x,y,z) coordinates of 4th vertex of the tetrahedron
//! \return Determinant of the Jacobian of transformation of physical
//!   tetrahedron to reference (xi, eta, zeta) space
// *****************************************************************************
{
  std::array< real, 3 > ba{{ v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] }},
                        ca{{ v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2] }},
                        da{{ v4[0]-v1[0], v4[1]-v1[1], v4[2]-v1[2] }};
  return triple( ba, ca, da );
}

std::array< std::array< tk::real, 3 >, 3 >
tk::inverseJacobian( const std::array< real, 3 >& v1,
                     const std::array< real, 3 >& v2,
                     const std::array< real, 3 >& v3,
                     const std::array< real, 3 >& v4 )
// *****************************************************************************
// Compute the inverse of the Jacobian of a coordinate transformation over a
// tetrahedron
//! \param[in] v1 (x,y,z) coordinates of 1st vertex of the tetrahedron
//! \param[in] v2 (x,y,z) coordinates of 2nd vertex of the tetrahedron
//! \param[in] v3 (x,y,z) coordinates of 3rd vertex of the tetrahedron
//! \param[in] v4 (x,y,z) coordinates of 4th vertex of the tetrahedron
//! \return Inverse of the Jacobian of transformation of physical
//!   tetrahedron to reference (xi, eta, zeta) space
// *****************************************************************************
{
  std::array< std::array< real, 3 >, 3 > jacInv;

  auto detJ = Jacobian( v1, v2, v3, v4 );

  jacInv[0][0] =  (  (v3[1]-v1[1])*(v4[2]-v1[2])
                   - (v4[1]-v1[1])*(v3[2]-v1[2])) / detJ;
  jacInv[1][0] = -(  (v2[1]-v1[1])*(v4[2]-v1[2])
                   - (v4[1]-v1[1])*(v2[2]-v1[2])) / detJ;
  jacInv[2][0] =  (  (v2[1]-v1[1])*(v3[2]-v1[2])
                   - (v3[1]-v1[1])*(v2[2]-v1[2])) / detJ;

  jacInv[0][1] = -(  (v3[0]-v1[0])*(v4[2]-v1[2])
                   - (v4[0]-v1[0])*(v3[2]-v1[2])) / detJ;
  jacInv[1][1] =  (  (v2[0]-v1[0])*(v4[2]-v1[2])
                   - (v4[0]-v1[0])*(v2[2]-v1[2])) / detJ;
  jacInv[2][1] = -(  (v2[0]-v1[0])*(v3[2]-v1[2])
                   - (v3[0]-v1[0])*(v2[2]-v1[2])) / detJ;

  jacInv[0][2] =  (  (v3[0]-v1[0])*(v4[1]-v1[1])
                   - (v4[0]-v1[0])*(v3[1]-v1[1])) / detJ;
  jacInv[1][2] = -(  (v2[0]-v1[0])*(v4[1]-v1[1])
                   - (v4[0]-v1[0])*(v2[1]-v1[1])) / detJ;
  jacInv[2][2] =  (  (v2[0]-v1[0])*(v3[1]-v1[1])
                   - (v3[0]-v1[0])*(v2[1]-v1[1])) / detJ;

  return jacInv;
}

tk::real
tk::determinant( const std::array< std::array< tk::real, 3 >, 3 >& a )
// *****************************************************************************
//  Compute the determinant of 3x3 matrix
//!  \param[in] a 3x3 matrix
//!  \return Determinant of the 3x3 matrix
// *****************************************************************************
{
  return ( a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]) );
}

std::array < tk::real, 3 >
tk::cramer( const std::array< std::array< tk::real, 3 >, 3>& a,
            const std::array< tk::real, 3 >& b )
// *****************************************************************************
//  Solve a 3x3 system of equations using Cramer's rule
//!  \param[in] a 3x3 lhs matrix
//!  \param[in] b 3x1 rhs matrix
//!  \return Array of solutions of the 3x3 system
// *****************************************************************************
{
  auto de = determinant( a );

  auto nu(0.0);
  std::array < real, 3 > x;

  nu = determinant( {{{{b[0], a[0][1], a[0][2]}},
                      {{b[1], a[1][1], a[1][2]}},
                      {{b[2], a[2][1], a[2][2]}}}} );
  x[0] = nu/de;

  nu = determinant( {{{{a[0][0], b[0], a[0][2]}},
                      {{a[1][0], b[1], a[1][2]}},
                      {{a[2][0], b[2], a[2][2]}}}} );
  x[1] = nu/de;

  nu = determinant( {{{{a[0][0], a[0][1], b[0]}},
                      {{a[1][0], a[1][1], b[1]}},
                      {{a[2][0], a[2][1], b[2]}}}} );
  x[2] = nu/de;

  return x;
}
