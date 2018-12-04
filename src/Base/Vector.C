// *****************************************************************************
/*!
  \file      src/Base/Vector.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************

#include <array>

#include "Vector.h"

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
              const std::array< real, 3 >& v2, real j )
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

tk::real
tk::Jacobian( const std::array< tk::real, 3 >& v1,
              const std::array< tk::real, 3 >& v2,
              const std::array< tk::real, 3 >& v3,
              const std::array< tk::real, 3 >& v4 )
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
  std::array< tk::real, 3 > ba{{ v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2] }},
                            ca{{ v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2] }},
                            da{{ v4[0]-v1[0], v4[1]-v1[1], v4[2]-v1[2] }};
  return tk::triple( ba, ca, da );
}

std::array< std::array< tk::real, 3 >, 3 >
tk::inverseJacobian( const std::array< tk::real, 3 >& v1,
                     const std::array< tk::real, 3 >& v2,
                     const std::array< tk::real, 3 >& v3,
                     const std::array< tk::real, 3 >& v4 )
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
  std::array< std::array< tk::real, 3 >, 3 > jacInv;

  auto detJ = tk::Jacobian( v1, v2, v3, v4 );

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
