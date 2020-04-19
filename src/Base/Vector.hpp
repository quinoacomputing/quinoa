// *****************************************************************************
/*!
  \file      src/Base/Vector.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2020 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************
#ifndef Vector_h
#define Vector_h

#include <array>

#include "Types.hpp"

namespace tk {

//! Flip sign of vector components
void
flip( std::array< real, 3 >& v ) noexcept;

//! Compute the cross-product of two vectors
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[out] rx x coordinate of the product vector
//! \param[out] ry y coordinate of the product vector
//! \param[out] rz z coordinate of the product vector
#pragma omp declare simd
inline void
cross( real v1x, real v1y, real v1z,
       real v2x, real v2y, real v2z,
       real& rx, real& ry, real& rz )
{
  rx = v1y*v2z - v2y*v1z;
  ry = v1z*v2x - v2z*v1x;
  rz = v1x*v2y - v2x*v1y;
}

//! Compute the cross-product of two vectors
std::array< real, 3 >
cross( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 );

//! Compute the cross-product of two vectors divided by a scalar
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[in] j The scalar to divide the product with
//! \param[out] rx x coordinate of the product vector
//! \param[out] ry y coordinate of the product vector
//! \param[out] rz z coordinate of the product vector
#pragma omp declare simd uniform(j)
inline void
crossdiv( real v1x, real v1y, real v1z,
          real v2x, real v2y, real v2z,
          real j,
          real& rx, real& ry, real& rz )
{
  cross( v1x, v1y, v1z, v2x, v2y, v2z, rx, ry, rz );
  rx /= j;
  ry /= j;
  rz /= j;
}

//! Compute the cross-product of two vectors divided by a scalar
std::array< real, 3 >
crossdiv( const std::array< real, 3 >& v1,
          const std::array< real, 3 >& v2,
          real j );

//! Compute the dot-product of two vectors
real
dot( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 )
noexcept;

//! Compute length of a vector
real
length( const std::array< real, 3 >& v ) noexcept;

//! Scale vector to unit length
void
unit( std::array< real, 3 >& v );

//! Compute the triple-product of three vectors
//! \param[in] v1x x coordinate of the 1st vector
//! \param[in] v1y y coordinate of the 1st vector
//! \param[in] v1z z coordinate of the 1st vector
//! \param[in] v2x x coordinate of the 2nd vector
//! \param[in] v2y y coordinate of the 2nd vector
//! \param[in] v2z z coordinate of the 2nd vector
//! \param[in] v3x x coordinate of the 3rd vector
//! \param[in] v3y y coordinate of the 3rd vector
//! \param[in] v3z z coordinate of the 3rd vector
//! \return Scalar value of the triple product
#pragma omp declare simd
inline tk::real
triple( real v1x, real v1y, real v1z,
        real v2x, real v2y, real v2z,
        real v3x, real v3y, real v3z )
{
  real cx, cy, cz;
  cross( v2x, v2y, v2z, v3x, v3y, v3z, cx, cy, cz );
  return v1x*cx + v1y*cy + v1z*cz;
}

//! Compute the triple-product of three vectors
real
triple( const std::array< real, 3 >& v1,
        const std::array< real, 3 >& v2,
        const std::array< real, 3 >& v3 );

//! Rotate vector about X axis
std::array< real, 3 >
rotateX( const std::array< real, 3 >& v, real angle );

//! Rotate vector about Y axis
std::array< real, 3 >
rotateY( const std::array< real, 3 >& v, real angle );

//! Rotate vector about Z axis
std::array< real, 3 >
rotateZ( const std::array< real, 3 >& v, real angle );

//! \brief Compute the determinant of the Jacobian of a coordinate
//!  transformation over a tetrahedron
real
Jacobian( const std::array< real, 3 >& v1,
          const std::array< real, 3 >& v2,
          const std::array< real, 3 >& v3,
          const std::array< real, 3 >& v4 );

//! \brief Compute the inverse of the Jacobian of a coordinate transformation
//!   over a tetrahedron
std::array< std::array< real, 3 >, 3 >
inverseJacobian( const std::array< real, 3 >& v1,
                 const std::array< real, 3 >& v2,
                 const std::array< real, 3 >& v3,
                 const std::array< real, 3 >& v4 );

//! Compute the determinant of 3x3 matrix
tk::real
determinant( const std::array< std::array< tk::real, 3 >, 3 >& a );

//! Solve a 3x3 system of equations using Cramer's rule
std::array < tk::real, 3 >
cramer( const std::array< std::array< tk::real, 3 >, 3>& a,
        const std::array< tk::real, 3 >& b );

} // tk::

#endif // Vector_h
