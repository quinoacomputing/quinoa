// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace
{
   int passed = 0;
   int failed = 0;

   using namespace pegtl;

   template< typename Rule >
   bool test_parse( const std::string & string )
   {
      string_input< ascii_location > in( string );
      try {
	 dummy_parse< Rule >( in );
	 return true;
      }
      catch ( const parse_error & ) {
	 return false;
      }
   }

   template< typename Rule >
   void s( const std::string & i )
   {
      if ( ! test_parse< Rule >( i ) ) {
	 PEGTL_PRINT( __PRETTY_FUNCTION__ << " failed on input \"" << i << "\"" );
	 ++failed;
	 return;
      }
      ++passed;
   }

   template< typename Rule >
   void f( const std::string & i )
   {
      if ( test_parse< Rule >( i ) ) {
	 PEGTL_PRINT( __PRETTY_FUNCTION__ << " succeeded on input \"" << i << "\"" );
	 ++failed;
	 return;
      }
      ++passed;
   }

   template< typename Rule >
   void b( const std::string & i, const bool t )
   {
      return t ? s< Rule >( i ) : f< Rule >( i );
   }

   void simple_tests()
   {
      PEGTL_PRINT( __FUNCTION__ );
      {
         const std::string e;
         assert( e.empty() );
         s< eol >( e );
         s< success >( e );
         f< failure >( e );
      }
      {
         const std::string e( "abc" );
         s< at< one< 'a' > > >( e );
         f< at< one< 'b' > > >( e );
         f< not_at< one< 'a' > > >( e );
         s< not_at< one< 'b' > > >( e );
      }
   }

   template< char A >
   void atomic_tests()
   {
      PEGTL_PRINT( __FUNCTION__ << " " << int( A ) );

      for ( int i = 0; ! ( i & 0x80 ); ++i ) {
	 const char c = char( i );
	 const std::string t( 1, c );
	 assert( t.size() == 1 );
	 b< lf >( t, c == '\n' );
	 b< cr >( t, c == '\r' );
	 s< any >( t );
	 s< success >( t );
	 f< failure >( t );
	 s< opt< one< A > > >( t );
	 s< star< one< A > > >( t );
	 b< plus< one< A > > >( t, c == A );
	 //	 s< rep< 0, one< A > > >( t );
	 b< rep< 1, one< A > > >( t, c == A );
	 f< rep< 2, one< A > > >( t );
	 f< rep< 2, one< 'a' > > >( t );
	 f< rep< 2, one< 'Z' > > >( t );
	 f< rep< 2, one< '5' > > >( t );
	 f< rep< 2, one< '\0' > > >( t );
	 b< one< A > >( t, c == A );
	 b< at_one< A > >( t, c == A );
	 b< seq< at_one< A >, one< A > > >( t, c == A );
	 b< not_one< A > >( t, c != A );
	 b< at_not_one< A > >( t, c != A );
	 b< seq< at_not_one< A >, not_one< A > > >( t, c != A );
	 b< lower >( t, ::islower( c ) );
	 b< upper >( t, ::isupper( c ) );
	 b< alpha >( t, ::isalpha( c ) );
	 b< digit >( t, ::isdigit( c ) );
	 b< alnum >( t, ::isalnum( c ) );
	 b< blank >( t, ::isblank( c ) );
	 b< space >( t, ::isspace( c ) );
	 b< xdigit >( t, ::isxdigit( c ) );
	 b< string< A > >( t, c == A );
	 b< at_string< A > >( t, c == A );
	 b< seq< at_string< A >, string< A > > >( t, c == A );
	 b< one< A > >( t, c == A );
	 b< at_one< A > >( t, c == A );
	 b< seq< at_one< A >, one< A > > >( t, c == A );
	 b< range< A, A > >( t, c == A );
	 b< at_range< A, A > >( t, c == A );
	 b< seq< at_range< A, A >, range< A, A > > >( t, c == A );
	 b< not_one< A > >( t, c != A );
	 b< at_not_one< A > >( t, c != A );
	 b< seq< at_not_one< A >, not_one< A > > >( t, c != A );
	 b< not_range< A, A > >( t, c != A );
	 b< at_not_range< A, A > >( t, c != A );
	 b< seq< at_not_range< A, A >, not_range< A, A > > >( t, c != A );
      }
   }

   void combinator_tests()
   {
      const std::string t = "abcabc";

      s< range< 'a', 'b' > >( t );
      s< range< 'a', 'c' > >( t );
      f< range< 'A', 'B' > >( t );
      f< range< 'A', 'C' > >( t );

      typedef string< 'a', 'b', 'c' > abc;
      typedef string< 'a', 'b', 'd' > abd;
      typedef string< 'b', 'b', 'd' > bbd;

      s< abc >( t );
      s< seq< abc, abc, eof > >( t );
      s< seq< sor< bbd, abc >, abc, eof > >( t );
      s< seq< sor< abd, abc >, abc, eof > >( t );

      s< rep< 1, abc > >( t );
      s< rep< 2, abc > >( t );
      f< rep< 3, abc > >( t );
      s< star< abc > >( t );
      s< plus< abc > >( t );
      s< seq< plus< abc >, star< abc >, eol > >( t );
      s< seq< abc, plus< abc >, star< abc >, eol > >( t );
      s< seq< plus< abc >, eol > >( t );
      s< seq< plus< abc >, eof > >( t );
      s< seq< star< abc >, eol > >( t );
      s< seq< star< abc >, eof > >( t );

      const std::string d = "abcabcd";

      s< seq< plus< abc >, one< 'd' > > >( d );
      s< seq< star< abc >, one< 'd' > > >( d );
      s< seq< rep< 2, abc >, one< 'd' > > >( d );
      f< seq< rep< 2, abc >, abc > >( d );
      s< seq< abc, abc, one< 'd' > > >( d );

      s< must< abc > >( t );
      f< must< one< 'd' > > >( t );
      s< until< one< 'd' > > >( d );
      s< until< one< 'a' > > >( d );
      s< until< one< 'c' > > >( d );
      f< until< one< 'd' > > >( t );

      s< until< one< 'a' >, one< 'b' > > >( t );
      f< until< one< 'c' >, one< 'b' > > >( t );
      s< until< one< 'b' >, one< 'a' > > >( t );
      f< until< one< 'c' >, one< 'a' > > >( t );

      s< ifthen< one< 'a' >, one< 'b' > > >( t );
      f< ifthen< one< 'a' >, one< 'c' > > >( t );
      s< ifthen< one< 'b' >, one< 'z' > > >( t );

      s< seq< ifthen< one< 'a' >, one< 'b' > >, one< 'c' > > >( t );
      s< sor< ifthen< one< 'a' >, one< 'c' > >, one< 'a' > > >( t );
      s< seq< ifthen< one< 'b' >, one< 'z' > >, one< 'a' > > >( t );

      s< ifmust< one< 'a' >, one< 'b' > > >( t );
      f< ifmust< one< 'a' >, one< 'c' > > >( t );
      f< ifmust< one< 'b' >, one< 'z' > > >( t );

      s< seq< ifmust< one< 'a' >, one< 'b' > >, one< 'c' > > >( t );
      f< sor< ifmust< one< 'a' >, one< 'c' > >, one< 'a' > > >( t );
      s< sor< ifmust< one< 'b' >, one< 'z' > >, one< 'a' > > >( t );

      s< ifthenelse< one< 'a' >, one< 'b' >, one< 'z' > > >( t );
      f< ifthenelse< one< 'a' >, one< 'c' >, one< 'z' > > >( t );
      f< ifthenelse< one< 'b' >, one< 'z' >, one< 'z' > > >( t );
      s< ifthenelse< one< 'b' >, one< 'z' >, one< 'a' > > >( t );

      s< ifmustelse< one< 'a' >, one< 'b' >, one< 'z' > > >( t );
      f< ifmustelse< one< 'a' >, one< 'c' >, one< 'z' > > >( t );
      f< ifmustelse< one< 'b' >, one< 'z' >, one< 'z' > > >( t );
      s< ifmustelse< one< 'b' >, one< 'z' >, one< 'a' > > >( t );

      f< list< one< 'a' >, one< 'b' > > >( "" );
      f< list< one< 'a' >, one< 'b' > > >( "b" );
      s< list< one< 'a' >, one< 'b' > > >( "a" );
      s< list< one< 'a' >, one< 'b' > > >( "aba" );
      s< list< one< 'a' >, one< 'b' > > >( "ababa" );
      f< list< one< 'a' >, one< 'b' > > >( "ba" );
      f< list< one< 'a' >, one< 'b' > > >( " ba" );
      f< list< one< 'a' >, one< 'b' > > >( "ab " );
      f< list< one< 'a' >, one< 'b' > > >( " a" );
      f< list< one< 'a' >, one< 'b' > > >( " aba" );
   }

} //

int main()
{
   simple_tests();

   atomic_tests< 'a' >();
   atomic_tests< 'Z' >();
   atomic_tests< '5' >();
   atomic_tests< ' ' >();
   atomic_tests< '"' >();
   atomic_tests< '\n' >();
   atomic_tests< '\'' >();
   atomic_tests< '\0' >();

   combinator_tests();

   if ( failed ) {
      std::cout << "ERROR, tests passed=" << passed << " tests failed=" << failed << "!" << std::endl;
      return 1;
   }
   else {
      std::cout << "OK, all " << passed << " tests passed." << std::endl;
      return 0;
   }
}
