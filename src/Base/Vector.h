//******************************************************************************
/*!
  \file      src/Base/Vector.h
  \author    J. Bakosi
  \date      Tue 12 May 2015 10:41:09 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Vector algebra
  \details   Vector algebra.
*/
//******************************************************************************
#ifndef Vector_h
#define Vector_h

#include <array>

#include <Types.h>

namespace tk {

std::array< tk::real, 3 >
cross( const std::array< tk::real, 3 >& v1,
       const std::array< tk::real, 3 >& v2 )
//******************************************************************************
//! Compute the cross-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Cross-product
//! \author J. Bakosi
//******************************************************************************
{
  return {{ v1[1]*v2[2] - v2[1]*v1[2],
            v1[2]*v2[0] - v2[2]*v1[0],
            v1[0]*v2[1] - v2[0]*v1[1] }};
}

tk::real
dot( const std::array< tk::real, 3 >& v1,
     const std::array< tk::real, 3 >& v2 )
//******************************************************************************
//! Compute the dot-product of two vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \return Dot-product
//! \author J. Bakosi
//******************************************************************************
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

tk::real
triple( const std::array< tk::real, 3 >& v1,
        const std::array< tk::real, 3 >& v2,
        const std::array< tk::real, 3 >& v3 )
//******************************************************************************
//! Compute the triple-product of three vectors
//! \param[in] v1 1st vector
//! \param[in] v2 2nd vector
//! \param[in] v3 3rd vector
//! \return Triple-product
//! \author J. Bakosi
//******************************************************************************
{
  return dot( v1, cross(v2,v3) );
}

} // tk::

#endif // Vector_h
