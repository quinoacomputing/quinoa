// Copyright (c) 2008 Dr. Colin Hirsch 
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_PARSE_STRING_HH
#define COHI_PEGTL_PARSE_STRING_HH

namespace pegtl
{
   // Functions to parse input given as std::string.

   // Wrapper functions that add another convenience layer: instantiation
   // and initialisation of the input class. See file parse_generic.hh for
   // the wrapped functions.

   // The *_parse_string_* functions set up the parser to parse a string.

   template< typename TopRule, typename Location = dummy_location, typename ... States >
   void dummy_parse_string( const std::string & string, States && ... st )
   {
      string_input< Location > in( string );
      dummy_parse< TopRule >( in, std::forward< States >( st ) ... );
   }

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void basic_parse_string( const std::string & string, States && ... st )
   {
      string_input< Location > in( string );
      basic_parse< TopRule >( in, std::forward< States >( st ) ... );
   }

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void trace_parse_string( const bool trace, const std::string & string, States && ... st )
   {
      string_input< Location > in( string );
      trace_parse< TopRule >( trace, in, std::forward< States >( st ) ... );
   }

   // Please read the comment on the smart_parse_* functions in parse_generic.hh!

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void smart_parse_string( const bool trace, const std::string & string, States && ... st )
   {
      try {
	 dummy_parse_string< TopRule >( string, std::forward< States >( st ) ... );
      }
      catch ( const parse_error & ) {
	 trace_parse_string< TopRule >( trace, string, std::forward< States >( st ) ... );
      }
   }

} // pegtl

#endif
