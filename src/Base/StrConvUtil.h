//******************************************************************************
/*!
  \file      src/Base/StrConvUtil.h
  \author    J. Bakosi
  \date      Thu 24 Jul 2014 10:25:07 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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

//! Operator << for adding (concatenating) T to a std::basic_string for lvalues
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( std::basic_string< Ch, Tr >& lhs, const T& e ) {
  std::stringstream ss;
  ss << lhs << e;
  lhs = ss.str();
  return lhs;
}

//! Operator << for adding (concatenating) T to a std::basic_string for rvalues
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( std::basic_string< Ch, Tr >&& lhs, const T& e ) {
  return lhs << e;
}

} // tk::

#endif // StrConvUtil_h
