//******************************************************************************
/*!
  \file      src/Base/StrConvUtil.h
  \author    J. Bakosi
  \date      Sat 05 Apr 2014 01:29:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     String conversion utilities
  \details   String conversion utilities
*/
//******************************************************************************
#ifndef StrConvUtil_h
#define StrConvUtil_h

#include <sstream>

namespace tk {

//! Operator << for writing T (casting to unsigned int) to output streams
template< typename T, typename Ch, typename Tr >
std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& e ) {
  os << static_cast< unsigned int >( e );
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
