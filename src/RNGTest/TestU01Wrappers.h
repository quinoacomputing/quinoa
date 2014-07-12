//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Wrappers.h
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 09:04:00 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     TestU01 global-scope wrappers
  \details   TestU01 global-scope wrappers. For more info on why these functions
             are in global-scope, see Main/RNGTest.C, for more info on why they
             are static inline, see http://stackoverflow.com/a/9399539
*/
//******************************************************************************
#ifndef TestU01Wrappers_h
#define TestU01Wrappers_h

namespace rngtest {

extern std::map< tk::ctr::RawRNGType, tk::RNG > g_rng;

template< tk::ctr::RawRNGType id >
static inline double uniform( void*, void* )
//******************************************************************************
//  Global-scope TestU01 uniform RNG wrapper
//! \author  J. Bakosi
//******************************************************************************
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
//******************************************************************************
//  Global-scope TestU01 uniform RNG bits wrapper
//! \author  J. Bakosi
//******************************************************************************
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
//******************************************************************************
//  Global-scope create TestU01 external generator
//! \author  J. Bakosi
//******************************************************************************
{
  return unif01_CreateExternGen01( const_cast<char*>(name.c_str()),
                                   uniform< id >, uniform_bits< id > );
}

} // rngtest::

#endif // TestU01Wrappers_h
