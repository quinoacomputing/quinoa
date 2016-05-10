// *****************************************************************************
/*!
  \file      src/Base/CharmUtil.h
  \author    J. Bakosi
  \date      Tue 10 May 2016 01:57:55 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ utilities
  \details   Charm++ utilities
*/
// *****************************************************************************
#ifndef CharmUtil_h
#define CharmUtil_h

#include "NoWarning/ice_and.h"

namespace tk {

//! Type trait querying whether T is a strongly typed enum
template< typename T >
using is_enum_class = typename boost::type_traits::ice_and<
                        std::is_enum< T >::value,
                        !std::is_convertible< T, uint8_t >::value >;

} // tk::

#endif // CharmUtil_h
