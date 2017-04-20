// *****************************************************************************
/*!
  \file      src/RNGTest/TestU01Wrappers.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     TestU01 global-scope wrappers
  \details   TestU01 global-scope wrappers. For more info on why these functions
    are in global-scope, see Main/RNGTest.C, for more info on why they are
    static inline, see http://stackoverflow.com/a/9399539.
*/
// *****************************************************************************
#ifndef TestU01Wrappers_h
#define TestU01Wrappers_h

#include <map>

#include "RNG.h"

namespace rngtest {

extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

template< tk::ctr::RawRNGType id >
static inline double uniform( void*, void* )
// *****************************************************************************
//  Global-scope TestU01 uniform RNG wrapper
//! \details This functions is used as an external uniform random number
//!   generator (i.e., external to TestU01) that is called by the TestU01
//!   library, hence the required signature and hence the global scope. It is
//!   templated on a unique integer corresponding to the RNG type enum defined
//!   by tk::ctr::RNGType. Templating on the id enables the compiler to generate
//!   a different wrapper for a different RNG facilitating simultaneous calls to
//!   any or all wrappers as they are unique functions.
//! \return Random number generated as a double-precision floating point value
//! \author J. Bakosi
// *****************************************************************************
{
  double r = 0.0;
  const auto rng = g_rng.find( id );
  if (rng != end(g_rng))
    rng->second.uniform( CkMyPe(), 1, &r );
  else
    Throw( "RNG not found" );
  return r;
}

template< tk::ctr::RawRNGType id >
static inline unsigned long uniform_bits( void*, void* )
// *****************************************************************************
//  Global-scope TestU01 uniform RNG bits wrapper
//! \details This functions is used as an external uniform random number
//!   generator (i.e., external to TestU01) that is called by the TestU01
//!   library, hence the required signature and hence the global scope. It is
//!   templated on a unique integer corresponding to the RNG type enum defined
//!   by tk::ctr::RNGType. Templating on the id enables the compiler to generate
//!   a different wrapper for a different RNG facilitating simultaneous calls to
//!   any or all wrappers as they are unique functions.
//! \return Random number generated as a unsigned long integer value
//! \author J. Bakosi
// *****************************************************************************
{
  double r = 0.0;
  const auto rng = g_rng.find( id );
  if (rng != end(g_rng))
    rng->second.uniform( CkMyPe(), 1, &r );
  else
    Throw( "RNG not found" );
  return static_cast<unsigned long>(r * unif01_NORM32);
}

template< tk::ctr::RawRNGType id >
static inline unif01_Gen* createTestU01Gen( const std::string& name )
// *****************************************************************************
//  Global-scope create TestU01 external generator
//! \details This is the function that contains the TestU01 library call to
//!   register a TestU01-external random number generator that later can be
//!   subjected to the TestU01 batteries. It ties the unique global-scope
//!   wrappers templated on the unique RNG id, thus TestU01 will see them as
//!   different external generators.
//! \param[in] name Random number generator name
//! \author J. Bakosi
// *****************************************************************************
{
  return unif01_CreateExternGen01( const_cast<char*>(name.c_str()),
                                   uniform< id >, uniform_bits< id > );
}

} // rngtest::

#endif // TestU01Wrappers_h
