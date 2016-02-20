//******************************************************************************
/*!
  \file      src/Base/CharmUtil.h
  \author    J. Bakosi
  \date      Wed 31 Dec 2014 05:01:37 PM MST
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Charm++ utilities
  \details   Charm++ utilities
*/
//******************************************************************************
#ifndef CharmUtil_h
#define CharmUtil_h

#include <boost/type_traits/detail/ice_and.hpp>

namespace tk {

//! Type trait querying whether T is a strongly typed enum
template< typename T >
using is_enum_class = typename boost::type_traits::ice_and<
                        std::is_enum< T >::value,
                        !std::is_convertible< T, uint8_t >::value >;

} // tk::

#endif // CharmUtil_h
