//******************************************************************************
/*!
  \file      src/Control/CustomPEGTLCapture.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 09:26:48 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Custom capture struct for PEGTL with std::transform
  \details   Custom capture struct for PEGTL with std::transform
*/
//******************************************************************************
#ifndef CustomPEGTLCapture_h
#define CustomPEGTLCapture_h

#include <pegtl.hh>

namespace tk {
namespace pegtl {

//! Custom capture based on PEGTL's capture, that applies the function convert
//! (e.g., toupper, tolower) on the map's value at match, see also
//! pegtl/rules_action.hh
template< unsigned Key, int (*convert)(int) >
struct capture
{
   typedef capture key_type;

   template< typename Print >
   static void prepare( Print & st )
   { st.template update< capture >( '\\' + ::pegtl::to_string( Key ), true ); }

   template< typename Map >
   static void apply( const std::string & s, Map & m ) { m[ Key ] = s; }

   template< bool Must, typename Input, typename Debug, typename Map >
   static bool match( Input & in, Debug &, const Map & m ) {
     typename Input::template marker< Must > p( in );

     const typename Map::const_iterator i = m.find( Key );

     if ( i == m.end() ) return p( false );
     std::string s = i->second;
     std::transform( begin(s), end(s), begin(s), convert );

     for ( size_t j = 0; j < s.size(); ++j ) {
       ::pegtl::character< Input > c( in );
       if ( ! c( c == s[ j ] ) ) return p( false );
     }
     return p( true );
   }
};

} // pegtl::
} // tk::

#endif // CustomPEGTLCapture_h
