//******************************************************************************
/*!
  \file      src/Base/Constants.h
  \author    J. Bakosi
  \date      Thu 11 Dec 2014 07:54:11 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Constants
  \details   Constants.
*/
//******************************************************************************
#ifndef Constants_h
#define Constants_h

#include <boost/math/constants/constants.hpp>

namespace tk {

//! Return the value of Pi

template< typename T >
constexpr T pi() {
  return boost::math::constants::pi< T >();
}

} // tk::

#endif // Constants_h
