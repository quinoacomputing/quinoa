// *****************************************************************************
/*!
  \file      src/Base/StrConvUtil.h
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     String conversion utilities
  \details   Various string conversion utilities.
*/
// *****************************************************************************
#ifndef StrConvUtil_h
#define StrConvUtil_h

#include <sstream>
#include <type_traits>

namespace tk {

//! Delegate operator << to default for writing non-enums to output streams.
//! \param[in] os Output stream into which t is written
//! \param[in] e  Value of arbitrary non-enum-class type to write to stream
//! \return Updated output stream for chain-use of the operator
template< typename T, typename Ch, typename Tr >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& e ) {
  if constexpr( std::is_enum_v<T> )
    os << static_cast< unsigned int >( e );
  else
    os << e;
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

} // tk::

#endif // StrConvUtil_h
