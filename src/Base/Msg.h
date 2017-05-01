// *****************************************************************************
/*!
  \file      src/Base/Msg.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Custom Charm++ message types definitions.
  \details   This file defines some custom Charm++ message types. This is
    somewhat experimental, not well tested, and is not currently in use. The
    best is to avoid this kind of low-level messages as this can be error-prone.
*/
// *****************************************************************************
#ifndef Msg_h
#define Msg_h

#include <cstddef>
#include <cstring>
#include <ostream>
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

#include "Exception.h"
#include "msg.decl.h"

namespace tk {

static const int CSTYLE_STRLEN = 1024;

//! Operator << for writing std::vector< T > to output streams; for debugging.
//! \param[in] os Output stream to write to
//! \param[in] t std::vector of arbitrary type T to write to stream
//! \return Updated output stream
//! \author J. Bakosi
template< typename T, typename Ch, typename Tr >
inline std::basic_ostream< Ch, Tr >&
operator<< ( std::basic_ostream< Ch, Tr >& os, const std::vector< T >& t ) {
  for (const auto& v : t) os << "'" << v << "' s:" << v.size() << " ";
  return os;
}

//! Charm++ message type for sending a single T, T must be POD.
//! \author J. Bakosi
template< typename T >
struct Msg : public CMessage_Msg< T > {
  using value_type = T;
  //! Constructor
  explicit Msg( value_type&& v ) : value( std::move(v) ) {}
  //! Value accessor
  value_type get() const { return value; }
  //! Data
  value_type value;
};

//! Charm++ message type for sending a string of strings separated by ';'.
//! \author J. Bakosi
struct StringsMsg : public CMessage_StringsMsg {
  using value_type = std::vector< std::string >;
  //! Constructor
  explicit StringsMsg( const std::string& v )
  { strncpy( strings, v.c_str(), CSTYLE_STRLEN ); }
  //! Value accessor splitting the migrated C-style string at semi-colons
  value_type get() {
    std::string s( strings );
    boost::char_separator<char> sep(";");
    boost::tokenizer< boost::char_separator<char> > tokens( s, sep );
    value_type v;
    for (const auto& t : tokens) v.push_back( t ); 
    return v;
  }
  //! Data
  char strings[CSTYLE_STRLEN];
};

//! Charm++ message type for sending a vector of strings of strings separated by
//! semi-colons. The size of the vector must be determined at compile-time and
//! comes in via the template argument. This avoids using variable-size
//! messages, which is more complicated, more error-prone, and most likely less
//! amenable to optimizations. The constructor and the value accessor interface
//! with std::vector< std::string >, but converted to POD so that the Charm++
//! runtime system easily serialize the data and migrate it over the network.
//! \author J. Bakosi
template< std::size_t Size = 1 >
struct VecStrsMsg : public CMessage_VecStrsMsg< Size > {
  using value_type = std::vector< std::vector< std::string > >;
  //! Constructor
  explicit VecStrsMsg( const std::vector< std::string >& v ) {
    Assert( Size <= v.size(), "Size of input vector to VecStrsMsg is larger "
                              "than specified in template argument" );
    for (std::size_t i=0; i<Size; ++i)
      strncpy( strings[i], v[i].c_str(), CSTYLE_STRLEN );
  }
  //! Value accessor
  value_type get() {
    value_type v( Size );
    for (std::size_t i=0; i<Size; ++i) {
      std::string s( strings[i] );
      boost::char_separator<char> sep(";");
      boost::tokenizer< boost::char_separator<char> > tokens( s, sep );
      for (const auto& t : tokens) v[i].push_back( t ); 
    }
    return v;
  }
  //! Data
  char strings[Size][CSTYLE_STRLEN];
};

//! Wait for and return future. Concept: Msg must have public function get()
//! returning Msg::value_type.
//! \author J. Bakosi
template< typename Msg >
typename Msg::value_type waitfor( const CkFuture& f ) {
  Msg* m = static_cast< Msg* >( CkWaitFuture( f ) );
  typename Msg::value_type value( m->get() );
  delete m;
  return value;
}

} // tk::

#define CK_TEMPLATES_ONLY
#include "msg.def.h"
#undef CK_TEMPLATES_ONLY

#endif // Msg_h
