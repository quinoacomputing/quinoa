// *****************************************************************************
/*!
  \file      src/Base/PrintUtil.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     String conversion utilities
  \details   Various string conversion utilities.
*/
// *****************************************************************************
#ifndef PrintUtil_h
#define PrintUtil_h

#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "Has.hpp"

namespace tk {

//! Operator << for writing an enum class to an output stream
//! \param[in] os Output stream into to write to
//! \param[in] e Value of enum-class type to write to stream
//! \return Updated output stream for chain-use of the operator
template< typename Enum, typename Ch, typename Tr,
          typename std::enable_if_t< std::is_enum_v<Enum>, int > = 0 >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const Enum& e ) {
  os << static_cast< unsigned int >( e );
  return os;
}

//! Operator << for writing a std::vector to an output stream
//! \param[in] os Output stream to write to
//! \param[in] v Vector to write to stream
//! \return Updated output stream for chain-use of the operator
template< class T, typename Ch, typename Tr >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const std::vector< T >& v ) {
  os << std::boolalpha;
  os << "[ ";
  for (const auto& p : v) os << p << ' ';
  os << ']';
  return os;
}

//! Operator << for writing an std::map to an output stream
//! \param[in] os Output stream to write to
//! \param[in] m Map to write to stream
//! \return Updated output stream for chain-use of the operator
template< typename Ch, typename Tr,
          class Key, class Value, class Compare = std::less< Key > >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os,
             const std::map< Key, Value, Compare >& m )
{
  if constexpr( tk::HasTypedef_i_am_tagged_tuple_v< Value > )
    for (const auto& [k,v] : m) os << '(' << k << ") : { " << v << "} ";
  else
    for (const auto& [k,v] : m) os << '(' << k << ") " << v << ' ';
  return os;
}

//! Operator << for adding (concatenating) T to a std::basic_string for lvalues.
//! \param[in] lhs Output std::basic_string into which e is written
//! \param[in] e Value of arbitrary type to write to string
//! \return Updated string
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( std::basic_string< Ch, Tr >& lhs, const T& e ) {
  std::stringstream ss;
  ss << lhs << e;
  lhs = ss.str();
  return lhs;
}

//! Operator << for adding (concatenating) T to a std::basic_string for rvalues.
//! \param[in] lhs Output std::basic_string into which e is written
//! \param[in] e Value of arbitrary type to write to string
//! \return Updated string
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( std::basic_string< Ch, Tr >&& lhs, const T& e ) {
  return lhs << e;
}

//!  Clean up whitespaces and format a long string into multiple lines
std::string
splitLines( std::string str,
            std::string indent,
            const std::string& name = "",
            std::size_t width = 80 );

// Calculate base log file name
std::string
baselogname( const std::string& executable );

//! Construct log file name
std::string
logname( const std::string& executable, int numrestart = 0 );

} // tk::

#endif // PrintUtil_h
