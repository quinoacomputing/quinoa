// *****************************************************************************
/*!
  \file      src/Base/Vector.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Vector algebra
  \details   Vector algebra.
*/
// *****************************************************************************
#ifndef Vector_h
#define Vector_h

#include <array>

namespace tk {

template< typename T >
std::array< T, 3 >
cross( const std::array< T, 3 >& v1, const std::array< T, 3 >& v2 )
// *****************************************************************************
//! Compute the cross-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Cross-product
//! \author J. Bakosi
// *****************************************************************************
{
  return {{ v1[1]*v2[2] - v2[1]*v1[2],
            v1[2]*v2[0] - v2[2]*v1[0],
            v1[0]*v2[1] - v2[0]*v1[1] }};
}

template< typename T >
std::array< T, 3 >
crossdiv( const std::array< T, 3 >& v1, const std::array< T, 3 >& v2, T j )
// *****************************************************************************
//! Compute the cross-product of two vectors divided by a scalar
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] j Scalar to divide each component by
//! \return Cross-product divided by scalar
//! \author J. Bakosi
// *****************************************************************************
{
  return {{ (v1[1]*v2[2] - v2[1]*v1[2]) / j,
            (v1[2]*v2[0] - v2[2]*v1[0]) / j,
            (v1[0]*v2[1] - v2[0]*v1[1]) / j }};
}

template< typename T >
T
dot( const std::array< T, 3 >& v1, const std::array< T, 3 >& v2 )
// *****************************************************************************
//! Compute the dot-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Dot-product
//! \author J. Bakosi
// *****************************************************************************
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

template< typename T >
T
triple( const std::array< T, 3 >& v1,
        const std::array< T, 3 >& v2,
        const std::array< T, 3 >& v3 )
// *****************************************************************************
//! Compute the triple-product of three vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] v3 3rd vector
//! \return Triple-product
//! \author J. Bakosi
// *****************************************************************************
{
  return dot( v1, cross(v2,v3) );
}

} // tk::

#endif // Vector_h
