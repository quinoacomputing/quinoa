//******************************************************************************
/*!
  \file      src/Base/CharmUtil.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 05:37:36 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Charm++ utilities
  \details   Charm++ utilities
*/
//******************************************************************************
#ifndef CharmUtil_h
#define CharmUtil_h

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

} // tk::

#endif // CharmUtil_h
