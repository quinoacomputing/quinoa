//******************************************************************************
/*!
  \file      src/Base/CharmUtil.h
  \author    J. Bakosi
  \date      Thu 03 Jul 2014 03:44:28 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Charm++ utilities
  \details   Charm++ utilities
*/
//******************************************************************************
#ifndef CharmUtil_h
#define CharmUtil_h

#include <boost/type_traits/detail/ice_and.hpp>

namespace tk {

//! Detect if T defines type "Proxy", credits go to Kerrek B at stackoverflow,
//! see http://stackoverflow.com/a/7235647
template< typename T >
struct HasProxy {
  private:
    typedef char                      yes;
    typedef struct { char array[2]; } no;
    template<typename C> static yes test(typename C::Proxy*);
    template<typename C> static no  test(...);
  public:
    static const bool value = sizeof(test<T>(0)) == sizeof(yes);
};

//! Type trait querying whether T is a strongly typed enum
template< typename T >
using is_enum_class = typename boost::type_traits::ice_and<
                        std::is_enum< T >::value,
                        !std::is_convertible< T, uint8_t >::value >;

} // tk::

#endif // CharmUtil_h
