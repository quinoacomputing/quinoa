// *****************************************************************************
/*!
  \file      src/UnitTest/TUTUtil.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Utilities for unit testing with the Template Unit Test library
  \details   Utilities for unit testing with the Template Unit Test library.
*/
// *****************************************************************************
#ifndef TUTUtil_h
#define TUTUtil_h

#include <limits>
#include <algorithm>

#include "NoWarning/tut.h"

namespace unittest {

//! \brief Ensure equality of all element of a vector of Ts (e.g., floating
//!   point numbers) up to some precision
//! \param[in] msg Message to output if the vectors are not equal
//! \param[in] a First vector to compare
//! \param[in] b Second vector to compare
//! \param[in] prec Optional precision
template< typename T >
void veceq( const std::string& msg,
            const std::vector< T >& a,
            const std::vector< T >& b,
            tk::real prec = std::numeric_limits< T >::epsilon() )
{
  std::equal( a.cbegin(), a.cend(), b.cbegin(),
              [ &msg, &prec ]( T s, T d )
              { tut::ensure_equals( msg, s, d, prec ); return true; } );
}

//! \brief Ensure equality of all element of a array of Ts (e.g., floating
//!   point numbers) up to some precision
//! \param[in] msg Message to output if the arrays are not equal
//! \param[in] a First array to compare
//! \param[in] b Second array to compare
//! \param[in] prec Optional precision
template< typename T, std::size_t N >
void veceq( const std::string& msg,
            const std::array< T, N >& a,
            const std::array< T, N >& b,
            tk::real prec = std::numeric_limits< T >::epsilon() )
{
  std::equal( a.cbegin(), a.cend(), b.cbegin(),
              [ &msg, &prec ]( T s, T d )
              { tut::ensure_equals( msg, s, d, prec ); return true; } );
}

} // unittest::

#endif // TUTUtil_h
