//******************************************************************************
/*!
  \file      src/RNG/TestU01Wrap.h
  \author    J. Bakosi
  \date      Mon 25 Nov 2013 11:08:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Interfacing the TestU01 random number generator test suite
  \details   Interfacing the TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01Wrap_h
#define TestU01Wrap_h

#include <Battery.h>

namespace rngtest {

//! Custom deleter binding a TestU01 pointer to its TestU01 deleter
template< class Ptr, void (*Del)(Ptr *) >
struct Eraser {
  void operator()( Ptr* gen ) {
    Del( gen );
  }
};

//! TestU01 pointer type with a custom deleter
template< class Ptr, void (*Del)(Ptr *) >
using TestU01Ptr = std::unique_ptr< Ptr, Eraser< Ptr, Del > >;

} // rngtest::

#endif // TestU01Wrap_h
