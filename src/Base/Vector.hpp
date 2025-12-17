// *****************************************************************************
/*!
  \file      src/Base/Vector.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************
#ifndef Vector_h
#define Vector_h

#include <array>
#include <cmath>
#include <vector>
#include <cblas.h>

#include "Types.hpp"
#include "Exception.hpp"

// ignore old-style-casts required for lapack/blas calls
#if defined(__clang__)
  #pragma clang diagnostic ignored "-Wold-style-cast"
#endif

// Lapacke forward declarations
extern "C" {

using lapack_int = long;

#define LAPACK_ROW_MAJOR 101
#define LAPACK_COL_MAJOR 102

extern lapack_int LAPACKE_dgetrf( int, lapack_int, lapack_int, double*,
  lapack_int, lapack_int* );
extern lapack_int LAPACKE_dgetri( int, lapack_int, double*, lapack_int,
  const lapack_int* );

}

namespace tk {

//! Flip sign of vector components
//! \param[in] v Vector whose components to multiply by -1.0
inline void
flip( std::array< real, 3 >& v )
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

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
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Cross-product
inline std::array< real, 3 >
cross( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 )
{
  real rx, ry, rz;
  cross( v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], rx, ry, rz );
  return { std::move(rx), std::move(ry), std::move(rz) };
}

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
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] j Scalar to divide each component by
//! \return Cross-product divided by scalar
inline std::array< real, 3 >
crossdiv( const std::array< real, 3 >& v1,
          const std::array< real, 3 >& v2,
          real j )
{
  return {{ (v1[1]*v2[2] - v2[1]*v1[2]) / j,
            (v1[2]*v2[0] - v2[2]*v1[0]) / j,
            (v1[0]*v2[1] - v2[0]*v1[1]) / j }};
}

//! Compute the dot-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Dot-product
inline real
dot( const std::array< real, 3 >& v1, const std::array< real, 3 >& v2 )
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

//! Compute the dot-product of a matrix and a vector
//! \param[in] m Matrix
//! \param[in] v vector
//! \return Dot-product
inline std::array< real, 3 >
matvec(
  const std::array< std::array< real, 3 >, 3 >& m,
  const std::array< real, 3 >& v )
{
  std::array< real, 3 > mv{0, 0, 0};
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      mv[i] += m[i][j]*v[j];
  }

  return mv;
}

//! Compute length of a vector
//! \param[in] x X coordinate of vector
//! \param[in] y Y coordinate of vector
//! \param[in] z Z coordinate of vector
//! \return length
#pragma omp declare simd
inline real
length( real x, real y, real z )
{
  return std::sqrt( x*x + y*y + z*z );
}

//! Compute length of a vector
//! \param[in] v vector
//! \return length
inline real
length( const std::array< real, 3 >& v )
{
  return std::sqrt( dot(v,v) );
}

//! Scale vector to unit length
//! \param[in,out] v Vector to normalize
inline void
unit( std::array< real, 3 >& v ) noexcept(ndebug)
{
  auto len = length( v );
  Assert( len > std::numeric_limits< tk::real >::epsilon(), "div by zero" );
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

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
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] v3 3rd vector
//! \return Triple-product
inline real
triple( const std::array< real, 3 >& v1,
        const std::array< real, 3 >& v2,
        const std::array< real, 3 >& v3 )
{
  return dot( v1, cross(v2,v3) );
}

//! Rotate vector about X axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateX( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ 1.0,         0.0,          0.0 }},
        {{ 0.0,   cos(angle), -sin(angle) }},
        {{ 0.0,   sin(angle),  cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! Rotate vector about Y axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateY( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle),  0.0, sin(angle) }},
        {{ 0.0,         1.0,        0.0 }},
        {{ -sin(angle), 0.0, cos(angle) }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! Rotate vector about Z axis
//! \param[in] v Vector to rotate
//! \param[in] angle Angle to use to rotate with
//! \return Rotated vector
inline std::array< real, 3 >
rotateZ( const std::array< real, 3 >& v, real angle )
{
  using std::cos;  using std::sin;

  std::array< std::array< real, 3 >, 3 >
    R{{ {{ cos(angle), -sin(angle), 0.0 }},
        {{ sin(angle),  cos(angle), 0.0 }},
        {{ 0.0,         0.0,        1.0 }} }};

  return {{ dot(R[0],v), dot(R[1],v), dot(R[2],v) }};
}

//! \brief Compute the determinant of the Jacobian of a coordinate
//!  transformation over a tetrahedron
//! \param[in] v1 (x,y,z) coordinates of 1st vertex of the tetrahedron
//! \param[in] v2 (x,y,z) coordinates of 2nd vertex of the tetrahedron
//! \param[in] v3 (x,y,z) coordinates of 3rd vertex of the tetrahedron
//! \param[in] v4 (x,y,z) coordinates of 4th vertex of the tetrahedron
//! \return Determinant of the Jacobian of transformation of physical
//!   tetrahedron to reference (xi, eta, zeta) space
inline real
Jacobian( const std::array< real, 3 >& v1,
          const std::array< real, 3 >& v2,
          const std::array< real, 3 >& v3,
          const std::array< real, 3 >& v4 )
{
  std::array< real, 3 > ba{{ v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] }},
                        ca{{ v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2] }},
                        da{{ v4[0]-v1[0], v4[1]-v1[1], v4[2]-v1[2] }};
  return triple( ba, ca, da );
}

//! \brief Compute the inverse of the Jacobian of a coordinate transformation
//!   over a tetrahedron
//! \param[in] v1 (x,y,z) coordinates of 1st vertex of the tetrahedron
//! \param[in] v2 (x,y,z) coordinates of 2nd vertex of the tetrahedron
//! \param[in] v3 (x,y,z) coordinates of 3rd vertex of the tetrahedron
//! \param[in] v4 (x,y,z) coordinates of 4th vertex of the tetrahedron
//! \return Inverse of the Jacobian of transformation of physical
//!   tetrahedron to reference (xi, eta, zeta) space
inline std::array< std::array< real, 3 >, 3 >
inverseJacobian( const std::array< real, 3 >& v1,
                 const std::array< real, 3 >& v2,
                 const std::array< real, 3 >& v3,
                 const std::array< real, 3 >& v4 )
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

//! Compute the determinant of 3x3 matrix
//!  \param[in] a 3x3 matrix
//!  \return Determinant of the 3x3 matrix
inline tk::real
determinant( const std::array< std::array< tk::real, 3 >, 3 >& a )
{
  return ( a[0][0] * (a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1] * (a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2] * (a[1][0]*a[2][1]-a[1][1]*a[2][0]) );
}

//! Compute the inverse of 3x3 matrix
//!  \param[in] m 3x3 matrix
//!  \return Inverse of the 3x3 matrix
inline std::array< std::array< tk::real, 3 >, 3 >
inverse( const std::array< std::array< tk::real, 3 >, 3 >& m )
{
  tk::real det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                 m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                 m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  tk::real invdet = 1.0 / det;

  std::array< std::array< tk::real, 3 >, 3 > minv;
  minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

  return minv;
}

//! Solve a 3x3 system of equations using Cramer's rule
//!  \param[in] a 3x3 lhs matrix
//!  \param[in] b 3x1 rhs matrix
//!  \return Array of solutions of the 3x3 system
inline std::array < tk::real, 3 >
cramer( const std::array< std::array< tk::real, 3 >, 3>& a,
        const std::array< tk::real, 3 >& b )
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

//! Move a point to a reference space given coordinates of origin of that space
//!  \param[in] origin Origin of reference frame to which point is to be moved
//!  \param[in,out] point Point that needs to be moved to reference frame
inline void
movePoint( const std::array< tk::real, 3 >& origin,
  std::array< tk::real, 3 >& point )
{
  for (std::size_t i=0; i<3; ++i)
    point[i] -= origin[i];
}

//! Calculate rotation matrix given three rotations in degrees
//!  \param[in] angles Angles in 3D space by which point is to be rotated
//!  \return Rotation matrix associated with rotations
inline std::array< std::array< tk::real, 3 >, 3 >
anglesToRotMat( const std::array< tk::real, 3 >& angles)
{
  // Convert angles to radian
  tk::real pi = 4.0*std::atan(1.0);
  auto a = angles[0] * pi/180.0;
  auto b = angles[1] * pi/180.0;
  auto c = angles[2] * pi/180.0;

  // Rotation matrix
  std::array< std::array< tk::real, 3 >, 3 > rotMat;
  rotMat[0][0] = cos(b)*cos(c);
  rotMat[0][1] = - cos(b)*sin(c);
  rotMat[0][2] = sin(b);

  rotMat[1][0] = sin(a)*sin(b)*cos(c) + cos(a)*sin(c);
  rotMat[1][1] = - sin(a)*sin(b)*sin(c) + cos(a)*cos(c);
  rotMat[1][2] = - sin(a)*cos(b);

  rotMat[2][0] = - cos(a)*sin(b)*cos(c) + sin(a)*sin(c);
  rotMat[2][1] = cos(a)*sin(b)*sin(c) + sin(a)*cos(c);
  rotMat[2][2] = cos(a)*cos(b);

  return rotMat;
}

//! Rotate a point in 3D space by specifying rotation angles in degrees
//!  \param[in] angles Angles in 3D space by which point is to be rotated
//!  \param[in,out] point Point that needs to be rotated
inline void
rotatePoint( const std::array< tk::real, 3 >& angles,
  std::array< tk::real, 3 >& point )
{
  // Rotation matrix
  auto rotMat = anglesToRotMat(angles);

  // Apply rotation
  std::array< tk::real, 3 > x{{0.0, 0.0, 0.0}};
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j) {
      x[i] += rotMat[i][j]*point[j];
    }
  }
  point = x;
}

//! \brief Get the Right Cauchy-Green strain tensor from the inverse deformation
//! gradient tensor.
//! \param[in] g Inverse deformation gradient tensor
//! \return Right Cauchy-Green tensor
inline std::array< std::array< real, 3 >, 3 >
getRightCauchyGreen(const std::array< std::array< real, 3 >, 3 >& g)
{
  // allocate matrices
  double G[9], C[9];

  // initialize c-matrices
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      G[i*3+j] = g[i][j];
  }

  // get g.g^T
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
    3, 3, 3, 1.0, G, 3, G, 3, 0.0, C, 3);

  // get inv(g.g^T)
  lapack_int ipiv[3];

  #ifndef NDEBUG
  lapack_int ierr =
  #endif
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, C, 3, ipiv);
  Assert(ierr==0, "Lapack error in LU factorization of g.g^T");

  #ifndef NDEBUG
  lapack_int jerr =
  #endif
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, C, 3, ipiv);
  Assert(jerr==0, "Lapack error in inverting g.g^T");

  // Output C as 2D array
  return {{ {C[0], C[1], C[2]},
            {C[3], C[4], C[5]},
            {C[6], C[7], C[8]} }};
}

//! \brief Get the Left Cauchy-Green strain tensor from the inverse deformation
//! gradient tensor.
//! \param[in] g Inverse deformation gradient tensor
//! \return Left Cauchy-Green tensor
inline std::array< std::array< real, 3 >, 3 >
getLeftCauchyGreen(const std::array< std::array< real, 3 >, 3 >& g)
{
  // allocate matrices
  double G[9], b[9];

  // initialize c-matrices
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      G[i*3+j] = g[i][j];
  }

  // get g^T.g
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    3, 3, 3, 1.0, G, 3, G, 3, 0.0, b, 3);

  // get inv(g^T.g)
  lapack_int ipiv[3];

  #ifndef NDEBUG
  lapack_int ierr =
  #endif
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 3, 3, b, 3, ipiv);
  Assert(ierr==0, "Lapack error in LU factorization of g^T.g");

  #ifndef NDEBUG
  lapack_int jerr =
  #endif
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, 3, b, 3, ipiv);
  Assert(jerr==0, "Lapack error in inverting g^T.g");

  // Output b as 2D array
  return {{ {b[0], b[1], b[2]},
            {b[3], b[4], b[5]},
            {b[6], b[7], b[8]} }};
}

//! \brief Get the volume-preserving part of the right Cauchy-Green strain
//!   tensor from the inverse deformation gradient tensor.
//! \param[in] g Inverse deformation gradient tensor
//! \return Volume-preserving part of the right Cauchy-Green strain tensor
inline std::array< std::array< tk::real, 3 >, 3 >
getIsochorRightCauchyGreen(const std::array< std::array< real, 3 >, 3 >& g)
{
  auto Ct = tk::getRightCauchyGreen(g);
  auto detC = std::pow(tk::determinant(Ct), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      Ct[i][j] /= detC;
  }

  return Ct;
}

//! \brief Get the deviatoric Hencky strain tensor from the inverse deformation
//! gradient tensor.
//! \param[in] g Inverse deformation gradient tensor
//! \return Deviatoric Hencky strain tensor
inline std::array< std::array< real, 3 >, 3 >
getDevHencky(const std::array< std::array< real, 3 >, 3 >& g)
{
  // Get right volm-preserving part of Cauchy-Green strain tensor
  auto C = getIsochorRightCauchyGreen(g);

  std::array< std::array< real, 3 >, 3 > devH{{{0,0,0}, {0,0,0}, {0,0,0}}};

  // Compute approximation of Hencky strain tensor from section 3.4 of
  // Barton, P. T. (2019). An interface-capturing Godunov method for the
  // simulation of compressible solid-fluid problems. Journal of Computational
  // Physics, 390, 25-50.

  // get g.g^T (i.e. inv(C))
  double G[9], CInv[9];
  for (std::size_t i=0; i<3; ++i) {
    for (std::size_t j=0; j<3; ++j)
      G[i*3+j] = g[i][j];
  }
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
    3, 3, 3, 1.0, G, 3, G, 3, 0.0, CInv, 3);

  // volume-preserving part
  auto detCInv = std::pow(tk::determinant(
    {{ { CInv[0], CInv[1], CInv[2] },
       { CInv[3], CInv[4], CInv[5] },
       { CInv[6], CInv[7], CInv[8] } }} ), 1.0/3.0);
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      CInv[3*i+j] /= detCInv;

  // Compute (C-CInv)/4
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      devH[i][j] = 0.25*(C[i][j]-CInv[3*i+j]);

  // Compute deviatoric part
  tk::real trH = devH[0][0] + devH[1][1] + devH[2][2];

  for (std::size_t i=0; i<3; ++i)
    devH[i][i] -= trH/3.0;

  // Output devH
  return devH;
}

//! \brief Get transpose of a 3 by 3 matrix
//! \param[in] mat matrix to be transposed
//! \return transposed matrix
inline std::array< std::array< tk::real, 3 >, 3 >
transpose3by3(std::array< std::array< tk::real, 3 >, 3 > mat)
{
  std::array< std::array< tk::real, 3 >, 3 > transMat;
  for (size_t i=0; i<3; ++i)
    for (size_t j=0; j<3; ++j)
      transMat[i][j] = mat[j][i];
  return transMat;
}

//! \brief Get rotation matrix in 2D array form given a normal
//! direction. The remaining directions are given by section 5.3.1 of
//! Miller, G. H., and P. Colella. "A conservative three-dimensional
//! Eulerian method for coupled solid–fluid shock capturing."
//! Journal of Computational Physics 183.1 (2002): 26-82.
//! \param[in] r Coordinates of the first basis vector r = (rx,ry,rz).
//! \return rotation matrix
inline std::array< std::array< tk::real, 3 >, 3 >
getRotationMatrix(const std::array< tk::real, 3 >& r)
{
  std::array< std::array< tk::real, 3 >, 3 > rotMat;
  tk::real rx = r[0];
  tk::real ry = r[1];
  tk::real rz = r[2];
  if (std::abs(ry+rz) <= std::abs(ry-rz))
  {
    tk::real norm = std::sqrt(2*(1-rx*ry-rx*rz-ry*rz));
    rotMat[0][0] = rx;
    rotMat[0][1] = ry;
    rotMat[0][2] = rz;
    rotMat[1][0] = (ry-rz)/norm;
    rotMat[1][1] = (rz-rx)/norm;
    rotMat[1][2] = (rx-ry)/norm;
    rotMat[2][0] = (rx*(ry+rz)-ry*ry-rz*rz)/norm;
    rotMat[2][1] = (ry*(rx+rz)-rx*rx-rz*rz)/norm;
    rotMat[2][2] = (rz*(rx+ry)-rx*rx-ry*ry)/norm;
  }
  else
  {
    tk::real norm = std::sqrt(2*(1+rz*(ry-rx)+rx*ry));
    rotMat[0][0] = rx;
    rotMat[0][1] = ry;
    rotMat[0][2] = rz;
    rotMat[1][0] = (ry+rz)/norm;
    rotMat[1][1] = (rz-rx)/norm;
    rotMat[1][2] = (-rx-ry)/norm;
    rotMat[2][0] = (rx*(rz-ry)-ry*ry-rz*rz)/norm;
    rotMat[2][1] = (ry*(rx+rz)+rx*rx+rz*rz)/norm;
    rotMat[2][2] = (rz*(rx-ry)-rx*rx-ry*ry)/norm;
  }
  return rotMat;
}

//! \brief Rotate a second order tensor (e.g. a Strain/Stress matrix) from
//! the (x,y,z) to a new (r,s,t) coordinate system.
//! The directions are given by section 5.3.1 of
//! Miller, G. H., and P. Colella. "A conservative three-dimensional
//! Eulerian method for coupled solid–fluid shock capturing."
//! Journal of Computational Physics 183.1 (2002): 26-82.
//! \param[in] mat matrix to be rotated.
//! \param[in] r Coordinates of the first basis vector r = (rx,ry,rz).
//! \return rotated tensor
inline std::array< std::array< tk::real, 3 >, 3 >
rotateTensor(const std::array< std::array< tk::real, 3 >, 3 >& mat,
             const std::array< tk::real, 3 >& r )
{
  // define rotation matrix
  std::array< std::array< tk::real, 3 >, 3 > rotMatrix = getRotationMatrix(r);
  double rotMat[9];

  // define matrices
  double matAuxIn[9], matAuxOut[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
    {
      matAuxIn[i*3+j] = mat[i][j];
      rotMat[i*3+j] = rotMatrix[i][j];
    }

  // compute matAuxIn*rotMat and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
    3, 3, 3, 1.0, matAuxIn, 3, rotMat, 3, 0.0, matAuxOut, 3);

  // matAuxOut -> matAuxIn
  for (std::size_t i=0; i<9; i++)
  {
    matAuxIn[i]  = matAuxOut[i];
    matAuxOut[i] = 0.0;
  }

  // compute rotMat^T*matAuxIn and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, 3, 1.0, rotMat, 3, matAuxIn, 3, 0.0, matAuxOut, 3);

  // return matAuxOut as a 2D array
  return {{ {matAuxOut[0], matAuxOut[1], matAuxOut[2]},
            {matAuxOut[3], matAuxOut[4], matAuxOut[5]},
            {matAuxOut[6], matAuxOut[7], matAuxOut[8]} }};
}

//! \brief Rotate a second order tensor (e.g. a Strain/Stress matrix) from
//! the (x,y,z) to a new (r,s,t) coordinate system.
//! The directions are given by section 5.3.1 of
//! Miller, G. H., and P. Colella. "A conservative three-dimensional
//! Eulerian method for coupled solid–fluid shock capturing."
//! Journal of Computational Physics 183.1 (2002): 26-82.
//! \param[in] mat matrix to be rotated.
//! \param[in] r Coordinates of the first basis vector r = (rx,ry,rz).
//! \return rotated tensor
inline std::array< std::array< tk::real, 3 >, 3 >
unrotateTensor(const std::array< std::array< tk::real, 3 >, 3 >& mat,
               const std::array< tk::real, 3 >& r )
{
  // define rotation matrix
  std::array< std::array< tk::real, 3 >, 3 > rotMatrix = getRotationMatrix(r);
  double rotMat[9];

  // define matrices
  double matAuxIn[9], matAuxOut[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
    {
      matAuxIn[i*3+j] = mat[i][j];
      rotMat[i*3+j] = rotMatrix[i][j];
    }
  // compute matAuxIn*rotMat and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, 3, 1.0, matAuxIn, 3, rotMat, 3, 0.0, matAuxOut, 3);
  // matAuxOut -> matAuxIn
  for (std::size_t i=0; i<9; i++)
  {
    matAuxIn[i]  = matAuxOut[i];
    matAuxOut[i] = 0.0;
  }
  // compute rotMat^T*matAuxIn and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    3, 3, 3, 1.0, rotMat, 3, matAuxIn, 3, 0.0, matAuxOut, 3);
  // return matAuxOut as a 2D array
  return {{ {matAuxOut[0], matAuxOut[1], matAuxOut[2]},
            {matAuxOut[3], matAuxOut[4], matAuxOut[5]},
            {matAuxOut[6], matAuxOut[7], matAuxOut[8]} }};
}

//! \brief Rotate a vector (e.g. a velocity) from
//! the (x,y,z) to a new (r,s,t) coordinate system.
//! The directions are given by section 5.3.1 of
//! Miller, G. H., and P. Colella. "A conservative three-dimensional
//! Eulerian method for coupled solid–fluid shock capturing."
//! Journal of Computational Physics 183.1 (2002): 26-82.
//! \param[in] v Vector to be rotated.
//! \param[in] r Coordinates of the first basis vector r = (rx,ry,rz).
//! \return rotated vector
inline std::array< tk::real, 3 >
rotateVector( const std::array< tk::real, 3 >& v,
  const std::array< tk::real, 3 >& r )
{
  // define rotation matrix
  std::array< std::array< tk::real, 3 >, 3 > rotMat = getRotationMatrix(r);

  // return rotMat*v
  return matvec(rotMat,v);
}

//! \brief Unrotate a vector (e.g. a velocity) from
//! the (x,y,z) to a new (r,s,t) coordinate system.
//! The directions are given by section 5.3.1 of
//! Miller, G. H., and P. Colella. "A conservative three-dimensional
//! Eulerian method for coupled solid–fluid shock capturing."
//! Journal of Computational Physics 183.1 (2002): 26-82.
//! \param[in] v Vector to be rotated.
//! \param[in] r Coordinates of the first basis vector r = (rx,ry,rz).
//! \return rotated vector
inline std::array< tk::real, 3 >
unrotateVector( const std::array< tk::real, 3 >& v,
  const std::array< tk::real, 3 >& r )
{
  // define rotation matrix
  std::array< std::array< tk::real, 3 >, 3 > rotMat = getRotationMatrix(r);

  // return rotMat^T*v
  return matvec(transpose3by3(rotMat),v);
}

//! \brief Reflect a second order tensor (e.g. a Strain/Stress matrix)
//! \param[in] mat matrix to be rotated.
//! \param[in] reflectMat Reflection matrix
//! \return reflected tensor
inline std::array< std::array< tk::real, 3 >, 3 >
reflectTensor(const std::array< std::array< tk::real, 3 >, 3 >& mat,
              const std::array< std::array< tk::real, 3 >, 3 >& reflectMat)
{
  // define reflection matrix
  double refMat[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      refMat[i*3+j] = reflectMat[i][j];

  // define matAux (I need matrices as row major 1D arrays)
  double matAuxIn[9], matAuxOut[9];
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      matAuxIn[i*3+j] = mat[i][j];

  // compute matAuxIn*refMat and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
    3, 3, 3, 1.0, matAuxIn, 3, refMat, 3, 0.0, matAuxOut, 3);

  // matAuxOut -> matAuxIn
  for (std::size_t i=0; i<9; i++)
  {
    matAuxIn[i]  = matAuxOut[i];
    matAuxOut[i] = 0.0;
  }

  // compute refMat^T*matAuxIn and store it into matAuxOut
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
    3, 3, 3, 1.0, refMat, 3, matAuxIn, 3, 0.0, matAuxOut, 3);

  // return matAuxOut as a 2D array
  return {{ {matAuxOut[0], matAuxOut[1], matAuxOut[2]},
            {matAuxOut[3], matAuxOut[4], matAuxOut[5]},
            {matAuxOut[6], matAuxOut[7], matAuxOut[8]} }};
}

//! \brief Wrapper to multiply two 3x3 matrices
//! \param[in] A matrix 1
//! \param[in] B matrix 2
//! \return A*B
inline std::array< std::array< tk::real, 3 >, 3 >
matmult33(const std::array< std::array< tk::real, 3 >, 3 >& A,
          const std::array< std::array< tk::real, 3 >, 3 >& B)
{
  // Unrolled as cblas was giving issues
  auto AB = A;
  AB[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
  AB[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
  AB[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
  AB[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
  AB[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
  AB[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
  AB[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
  AB[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
  AB[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
  return AB;
}

//! \brief Multiply 2 quaternions
//! \param[in] a first quaternion
//! \param[in] b second quaternion
//! \return multiplied quaternion
inline std::array< tk::real, 4 >
quaternion_mult(const std::array< tk::real, 4 >& a,
                const std::array< tk::real, 4 >& b)
{
  std::array< tk::real, 3 > av{ a[1], a[2], a[3] };
  std::array< tk::real, 3 > bv{ b[1], b[2], b[3] };
  auto as = a[0];
  auto bs = b[0];
  auto abs = as*bs - dot(av, bv);
  auto abv = cross(av, bv);
  for (std::size_t i = 0; i < 3; ++i)
      abv[i] += as*bv[i] + bs*av[i];
  std::array< tk::real, 4 > ab{abs, abv[0], abv[1], abv[2] };
  return ab;
}

//! \brief Obtain the magnitude of a quaternion
//! \param[in] q quaternion
//! \return quaternion magnitude
inline tk::real
quaternion_mag(const std::array< tk::real, 4 >& q)
{
  std::array< tk::real, 3 > v{ q[1], q[2], q[3] };
  return std::sqrt(q[0]*q[0] + dot(v, v));
}

//! \brief Convert a rotation quaternion to a rotation matrix
//! \param[in] q quaternion
//! \return Rotation matrix
inline std::array< std::array< tk::real, 3 >, 3 >
qtoR(const std::array< tk::real, 4 >& q)
{
  std::array< std::array< tk::real, 3 >, 3 > R;
  std::array< tk::real, 3 > v{ q[1], q[2], q[3] };
  auto s = q[0];
  // The conversion for quaternion [s, v_i] is:
  // vx_{ij} = -e_{ijk} v_k
  // R_{ij} = delta_{ij} + 2*s*vx_{ij} + 2*vx_{ik}vx_{kj}
  // e_{ijk} is the Levi-Civita tensor
  std::array< std::array< tk::real, 3 >, 3 >  vx{
    {
      {0.0, -v[2], v[1]},
      {v[2], 0.0, -v[0]},
      {-v[1], v[0], 0.0},
    }};

  auto vxvx = matmult33(vx, vx);
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j){
      R[i][j] = 2.0*s*vx[i][j] + 2.0*vxvx[i][j];
      if (i == j) R[i][j] += 1.0;
    }
  return R;
}

//! \brief Convert a rotation matrix to a rotation quaternion
//! \param[in] R rotation amtrix
//! \return rotation quaternion
inline std::array< tk::real, 4 >
Rtoq(const std::array< std::array< tk::real, 3 >, 3 >& R)
{
  std::array< tk::real, 4 > q;

  auto tr = R[0][0] + R[1][1] + R[2][2];
  tk::real S, qw, qx, qy, qz;

  if (tr > 0) {
    S = sqrt(tr + 1.0) * 2; // S=4*qw
    qw = 0.25 * S;
    qx = (R[2][1] - R[1][2]) / S;
    qy = (R[0][2] - R[2][0]) / S;
    qz = (R[1][0] - R[0][1]) / S;
  } else if ((R[0][0] > R[1][1]) & (R[0][0] > R[2][2])) {
    S = sqrt(1.0 + R[0][0] - R[1][1] - R[2][2]) * 2;
    qw = (R[2][1] - R[1][2]) / S;
    qx = 0.25 * S;
    qy = (R[0][1] + R[1][0]) / S;
    qz = (R[0][2] + R[2][0]) / S;
  } else if (R[1][1] > R[2][2]) {
    S = sqrt(1.0 + R[1][1] - R[0][0] - R[2][2]) * 2;
    qw = (R[0][2] - R[2][0]) / S;
    qx = (R[0][1] + R[1][0]) / S;
    qy = 0.25 * S;
    qz = (R[1][2] + R[2][1]) / S;
  } else {
    S = sqrt(1.0 + R[2][2] - R[0][0] - R[1][1]) * 2;
    qw = (R[1][0] - R[0][1]) / S;
    qx = (R[0][2] + R[2][0]) / S;
    qy = (R[1][2] + R[2][1]) / S;
    qz = 0.25 * S;
  }

  q[0] = qw; q[1] = qx; q[2] = qy; q[3] = qz;
  return q;
}

//! Branchless Cholesky factorization of a 3x3 (SPD) matrix
//! \param[in] A Matrix to be factorized
//! \param[in,out] L Cholesky factorization
inline void chol3x3( const std::array< std::array< tk::real, 3 >, 3 >& A,
                     std::array< std::array< tk::real, 3 >, 3 >& L )
{
  // A is symmetric: [a b c; b d e; c e f]
  const tk::real a = A[0][0], b = A[0][1], c = A[0][2];
  const tk::real d = A[1][1], e = A[1][2], f = A[2][2];

  // Regularize (helps near-degenerate stencils)
  const tk::real tr = a + d + f;
  const tk::real eps = std::max<tk::real>(1e-30, tr * 1e-14);

  const tk::real aR = a + eps;
  const tk::real dR = d + eps;
  const tk::real fR = f + eps;

  const tk::real L00 = std::sqrt(aR);
  const tk::real L10 = b / L00;
  const tk::real L20 = c / L00;

  const tk::real t11 = dR - L10*L10;
  const tk::real L11 = std::sqrt( t11 > 0 ? t11 : eps );

  const tk::real L21 = (e - L10*L20) / L11;

  const tk::real t22 = fR - L20*L20 - L21*L21;
  const tk::real L22 = std::sqrt( t22 > 0 ? t22 : eps );

  L = {{ { L00, 0.0, 0.0 },
         { L10, L11, 0.0 },
         { L20, L21, L22 } }};
}

//! Solve 3x3 system using the Cholesky factorization the 3x3 (SPD) matrix
//! \param[in] L Cholesky factorization of 3x3 LHS matrix
//! \param[in] b RHS matrix
//! \return Solution of the 3x3 system Ax=b
inline std::array< tk::real,3 >
solve_chol3x3( const std::array< std::array< tk::real, 3 >, 3 >& L,
               const std::array< tk::real, 3 >& b )
{
  // Forward: L y = b
  const tk::real y0 = b[0] / L[0][0];
  const tk::real y1 = (b[1] - L[1][0]*y0) / L[1][1];
  const tk::real y2 = (b[2] - L[2][0]*y0 - L[2][1]*y1) / L[2][2];

  // Backward: L^T x = y
  tk::real x2 = y2 / L[2][2];
  tk::real x1 = (y1 - L[2][1]*x2) / L[1][1];
  tk::real x0 = (y0 - L[1][0]*x1 - L[2][0]*x2) / L[0][0];

  return {{ x0, x1, x2 }};
}

} // tk::

#endif // Vector_h
