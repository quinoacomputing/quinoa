//******************************************************************************
/*!
  \file      src/Base/StrConvUtil.h
  \author    J. Bakosi
  \date      Mon 02 Jun 2014 06:04:28 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     String conversion utilities
  \details   String conversion utilities
*/
//******************************************************************************
#ifndef StrConvUtil_h
#define StrConvUtil_h

#include <sstream>

namespace tk {

//! Operator << for writing enum class value to output streams
template< typename T, typename Ch, typename Tr,
          typename std::enable_if< std::is_enum<T>::value, int >::type = 0 >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& e ) {
  os << static_cast< unsigned int >( e );
  return os;
}

//! Delegate operator << to default for writing non-enums to output streams
template< typename T, typename Ch, typename Tr,
          typename std::enable_if< !std::is_enum<T>::value, int >::type = 0 >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& t ) {
  os << t;
  return os;
}

//! Operator << for adding (concatenating) T to a std::basic_strin
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( const std::basic_string< Ch, Tr >& lhs, const T& e ) {
  std::stringstream ss;
  ss << lhs << e;
  return ss.str();
}

//! Operator + for adding (concatenating) T to a std::basic_string
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator+ ( const std::basic_string< Ch, Tr >& lhs, const T& e ) {
  return lhs << e;      // implement in terms of operator<<
}

} // tk::

#endif // StrConvUtil_h
