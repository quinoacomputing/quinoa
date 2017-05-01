// *****************************************************************************
/*!
  \file      src/Base/StrConvUtil.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     String conversion utilities
  \details   Various string conversion utilities.
*/
// *****************************************************************************
#ifndef StrConvUtil_h
#define StrConvUtil_h

#include <sstream>

namespace tk {

//! Operator << for writing enum class value to output streams.
//! \param[in] os Output stream into which e is written
//! \param[in] e  Value of enum class to write to stream
//! \return Updated output stream for chain-use of the operator
//! \author J. Bakosi
template< typename T, typename Ch, typename Tr,
          typename std::enable_if< std::is_enum<T>::value, int >::type = 0 >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& e ) {
  os << static_cast< unsigned int >( e );
  return os;
}

//! Delegate operator << to default for writing non-enums to output streams.
//! \param[in] os Output stream into which t is written
//! \param[in] t  Value of arbitrary non-enum-class type to write to stream
//! \return Updated output stream for chain-use of the operator
//! \author J. Bakosi
template< typename T, typename Ch, typename Tr,
          typename std::enable_if< !std::is_enum<T>::value, int >::type = 0 >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const T& t ) {
  os << t;
  return os;
}

//! Operator << for adding (concatenating) T to a std::basic_string for lvalues.
//! \param[in] lhs Output std::basic_string into which e is written
//! \param[in] e Value of arbitrary type to write to string
//! \return Updated string
//! \author J. Bakosi
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
//! \author J. Bakosi
template< typename T, typename Ch, typename Tr >
std::basic_string< Ch, Tr >
operator<< ( std::basic_string< Ch, Tr >&& lhs, const T& e ) {
  return lhs << e;
}

} // tk::

#endif // StrConvUtil_h
