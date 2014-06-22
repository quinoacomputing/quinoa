//******************************************************************************
/*!
  \file      src/Base/Msg.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 10:25:18 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Custom Charm++ message types
  \details   Custom Charm++ message types
*/
//******************************************************************************
#ifndef Msg_h
#define Msg_h

#include <boost/tokenizer.hpp>
#include <msg.decl.h>

namespace tk {

//! Operator << for writing std::vector< T > to output streams
template< typename T, typename Ch, typename Tr >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const std::vector< T >& t ) {
  for (const auto& v : t) os << "'" << v << "' s:" << v.size() << " ";
  return os;
}

//! Charm++ message type for sending a single T, T must be POD
template< typename T >
struct Msg : public CMessage_Msg< T > {
  using value_type = T;
  explicit Msg( value_type&& v ) : value( std::move(v) ) {}
  value_type get() const { return value; }
  value_type value;
};

//! Charm++ message type for sending a string of strings separated by ';'
struct StringsMsg : public CMessage_StringsMsg {
  using value_type = std::vector< std::string >;
  explicit StringsMsg( const std::string& v ) { strncpy(str, v.c_str(), 1024); }
  value_type get() {
    std::string s( str );
    boost::char_separator<char> sep(";");
    boost::tokenizer< boost::char_separator<char> > tokens( s, sep );
    value_type v;
    for (const auto& t : tokens) v.push_back( t ); 
    return v;
  }
  char str[1024];
};

//! Wait for and return future
template< typename Msg >
typename Msg::value_type waitfor( const CkFuture& f ) {
  Msg* m = static_cast< Msg* >( CkWaitFuture( f ) );
  typename Msg::value_type value( m->get() );
  delete m;
  return value;
}

} // tk::

#define CK_TEMPLATES_ONLY
#include <msg.def.h>
#undef CK_TEMPLATES_ONLY

#endif // Msg_h
