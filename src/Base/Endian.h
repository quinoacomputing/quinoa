// *****************************************************************************
/*!
  \file      src/Control/Endian.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Swap endianness
  \details   Swap endianness. Thanks to http://stackoverflow.com/a/4956493 and
             http://stackoverflow.com/a/3522853.
*/
// *****************************************************************************
#ifndef Endian_h
#define Endian_h

#include <climits>

#include "Macro.h"

namespace tk {

//! Swap endianness of an integral type
//! \param[in] u Integral type to convert
//! \return Converted integer
template< typename T >
T swap_endian( T u ) {
  static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");

  union {
    T u;
    unsigned char u8[sizeof(T)];
  } source, dest;

  source.u = u;

  for (size_t k = 0; k < sizeof(T); k++)
    dest.u8[k] = source.u8[ sizeof(T) - k - 1 ];

  return dest.u;
}

#if defined(__clang__)
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wold-style-cast"
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

//! Swap endianness of a double
//! \param[in] u Double to convert
//! \return Converted double
template<>
double swap_endian( double u ) {
  uint64_t mem = swap_endian< uint64_t>( *(uint64_t*)&u );
  return *(double*)&mem;
}

#if defined(__clang__)
  #pragma clang diagnostic pop
#elif defined(STRICT_GNUC)
  #pragma GCC diagnostic pop
#endif

} // ::tk

#endif // Endian_h
