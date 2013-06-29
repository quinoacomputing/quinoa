// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_RULES_STRING_HH
#define COHI_PEGTL_RULES_STRING_HH

namespace pegtl
{
   struct any
   {
      typedef any key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< any >( ".", true );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 if ( in.eof() ) {
	    return false;
	 }
	 in.bump();
	 return true;
      }
   };

   template< int... Chars > struct one;

   template<>
   struct one<>
   {
      static bool i_match( const int )
      {
	 return false;
      }
   };

   template< int Char, int... Chars >
   struct one< Char, Chars... >
   {
      typedef one key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = ( sizeof ... ( Chars ) ) ? ( "\"[" + escaper< Char, Chars ... >::result() + "]\"" ) : ( '"' + escape( Char ) + '"' );
	 st.template update< one >( n, true );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 if ( in.eof() ) {
	    return false;
	 }
	 character< Input > h( in );
	 return h( i_match( h ) );
      }

      static bool i_match( const int i )
      {
	 return ( i == Char ) || one< Chars... >::i_match( i );
      }
   };

   template< int ... Chars >
   struct at_one
	 : at< one< Chars ... > > {};

   template< int ... Chars > struct not_one;

   template<>
   struct not_one<>
   { };

   template< int Char, int ... Chars >
   struct not_one< Char, Chars ... >
   {
      typedef not_one key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = std::string( "\"[^" ) + escaper< Char, Chars ... >::result() + "]\"";
	 st.template update< not_one >( n, true );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 if ( in.eof() ) {
	    return false;
	 }
	 character< Input > h( in );
	 return h( ! one< Char, Chars ... >::i_match( h ) );
      }
   };

   template< int ... Chars >
   struct at_not_one
	 : at< not_one< Chars ... > > {};

   template< int C, int D >
   struct range
   {
      typedef range key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = std::string( "\"[" ) + escape( C ) + "-" + escape( D ) + "]\"";
	 st.template update< range >( n, true );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 static_assert( C <= D, "pegtl: illegal expression range< C, D > where C is greater than D" );
	 if ( in.eof() ) {
	    return false;
	 }
	 character< Input > h( in );
	 return h( ( h >= C ) && ( h <= D ) );
      }
   };

   template< int C, int D >
   struct not_range
   {
      typedef not_range key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = std::string( "\"[^" ) + escape( C ) + "-" + escape( D ) + "]\"";
	 st.template update< not_range >( n, true );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 static_assert( C <= D, "pegtl: illegal expression not_range< C, D > where C is greater than D" );

	 if ( in.eof() ) {
	    return false;
	 }
	 character< Input > h( in );
	 return h( ( h < C ) || ( h > D ) );
      }
   };

   template< int C, int D >
   struct at_range
	 : at< range< C, D > > {};

   template< int C, int D >
   struct at_not_range
	 : at< not_range< C, D > > {};

   template< int... Chars > struct string;

   template<>
   struct string<>
	 : success
   {
      template< typename Print >
      static void i_insert( std::ostream &, Print & )
      { }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool i_match( Input &, Debug &, States && ... )
      {
	 return true;
      }
   };

   template< int Char, int ... Chars >
   struct string< Char, Chars ... >
   {
      typedef string key_type;

      template< typename Print >
      static void i_insert( std::ostream & o, Print & st )
      {
	 st.template insert< one< Char > >();
	 o << escape( Char );
	 string< Chars ... >::i_insert( o, st );
      }

      template< typename Print >
      static void prepare( Print & st )
      {
	 std::ostringstream o;
	 i_insert( o, st );
	 const std::string n = std::string( "\"" ) + o.str() + "\"";
	 st.template update< string >( n, true );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must > h( in );
	 return h( i_match< Must >( in, de, std::forward< States >( st ) ... ) );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool i_match( Input & in, Debug & de, States && ... st )
      {
	 return one< Char >::template match< Must >( in, de, std::forward< States >( st ) ... ) && pegtl::string< Chars ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... );
      }
   };

   template< int Char, int ... Chars >
   struct at_string
	 : at< string< Char, Chars ... > > {};

   struct lf
	 : one< '\n' > {};

   struct cr
	 : one< '\r' > {};

   struct crlf
	 : seq< cr, lf > {};

   struct eol
	 : sor< eof, crlf, lf, cr >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< eol >( "eol", true );
      }
   };

   struct digit
	 : range< '0', '9' >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< digit >( "digit", true );
      }
   };

   struct lower
	 : range< 'a', 'z' >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< lower >( "lower", true );
      }
   };

   struct upper
	 : range< 'A', 'Z' >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< upper >( "upper", true );
      }
   };

   struct alpha
	 : sor< lower, upper >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< alpha >( "alpha", true );
      }
   };

   struct alnum
	 : sor< alpha, digit >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< alpha >( "alnum", true );
      }
   };

   struct xdigit
	 : sor< digit, range< 'a', 'f' >, range< 'A', 'F' > >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< alpha >( "xdigit", true );
      }
   };

   struct blank
	 : one< ' ', '\t' >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< blank >( "blank", true );
      }
   };

   struct space
	 : one< ' ', '\n', '\r', '\t', '\v', '\f' >
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< space >( "space", true );
      }
   };

   typedef sor< one< '_' >, alpha > ident1;
   typedef sor< digit, ident1 > ident2;

   struct identifier
	 : seq< ident1, star< ident2 > > {};

   template< int Char, typename RulePadL, typename RulePadR = RulePadL >
   struct pad_one
	 : pad< one< Char >, RulePadL, RulePadR > {};

   struct space_until_eof
	 : until< eof, space > {};

   struct blank_until_eol
	 : until< eol, blank > {};

   struct shebang
	 : ifmust< string< '#', '!' >, until< eol > > {};

} // pegtl

#endif
