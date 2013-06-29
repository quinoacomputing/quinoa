// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <fstream>
#include <pegtl.hh>

namespace scheme
{
   // Nearly complete Scheme R6RS syntax check.
   // - Only ASCII input is handled correctly.
   // - The rule <u8> does not check the exact range.
   // - This is as faithful to the original as possible,
   //   not even obvious redundancies have been changed;
   //   some redundancies could be easily factored out.
   // - Sometimes, however, the order of rules was changed
   //   because while CFGs backtrack when a nondeterministic
   //   choice fails, PEGs are always deterministic because
   //   the order of rules implies a priority, and a choice
   //   that was once made is never changed later.
   // - Some additions were necessary in order to replace
   //   the otherwise separate tokenisation phase.

   using namespace pegtl::ascii;

   struct interlexeme_space;

   template< typename Lexeme >
   struct padded_lexeme
	 : pegtl::rule_base< padded_lexeme< Lexeme >, pegtl::pad< Lexeme, interlexeme_space > > {};

   template< int R > struct digit;

   struct hex_digit
	 : pegtl::xdigit {};

   template<>
   struct digit< 2 >
	 : pegtl::one< '0', '1' > {};

   template<>
   struct digit< 8 >
	 : pegtl::one< '0', '1', '2', '3', '4', '5', '6', '7' > {};

   template<>
   struct digit< 10 >
	 : pegtl::digit {};

   template<>
   struct digit< 16 >
	 : hex_digit {};

   template< int R > struct radix;

   template<>
   struct radix< 2 >
	 : pegtl::sor< pegtl::string< '#', b >,
		       pegtl::string< '#', B > > {};

   template<>
   struct radix< 8 >
	 : pegtl::sor< pegtl::string< '#', o >,
		       pegtl::string< '#', O > > {};

   template<>
   struct radix< 10 >
	 : pegtl::sor< pegtl::string< '#', d >,
		       pegtl::string< '#', D >,
		       pegtl::success > {};

   template<>
   struct radix< 16 >
	 : pegtl::sor< pegtl::string< '#', x >,
		       pegtl::string< '#', X > > {};

   struct exactness
	 : pegtl::sor< pegtl::string< '#', i >,
		       pegtl::string< '#', I >,
		       pegtl::string< '#', e >,
		       pegtl::string< '#', E >,
		       pegtl::success > {};

   struct sign
	 : pegtl::opt< pegtl::one< '+', '-' > > {};

   struct mantissa_width
	 : pegtl::opt< pegtl::seq< pegtl::one< '|' >, pegtl::plus< digit< 10 > > > > {};

   struct exponent_marker
	 : pegtl::one< e, E, s, S, f, F, d, D, l, L > {};

   struct suffix
	 : pegtl::opt< pegtl::seq< exponent_marker, sign, pegtl::plus< digit< 10 > > > > {};

   template< int R >
   struct prefix
	 : pegtl::sor< pegtl::seq< radix< R >, exactness >,
		       pegtl::seq< exactness, radix< R > > > {};

   template< int R >
   struct uinteger
	 : pegtl::plus< digit< R > > {};

   template< int R >
   struct decimal
	 : pegtl::failure {};

   template<>
   struct decimal< 10 >
	 : pegtl::sor< pegtl::seq< uinteger< 10 >, suffix >,
		       pegtl::seq< pegtl::one< '.' >, pegtl::plus< digit< 10 > >, suffix >,
		       pegtl::seq< pegtl::plus< digit< 10 > >, pegtl::one< '.' >, pegtl::star< digit< 10 > >, suffix >,
		       pegtl::seq< pegtl::plus< digit< 10 > >, pegtl::one< '.' >, suffix > > {};

   template< int R >
   struct ureal
	 : pegtl::sor< uinteger< R >,
		       pegtl::seq< uinteger< R >, pegtl::one< '/' >, uinteger< R > >,
		       pegtl::seq< decimal< R >, mantissa_width > > {};

   struct naninf
	 : pegtl::sor< pegtl::string< n, a, n, '.', '0' >,
		       pegtl::string< i, n, f, '.', '0' > > {};

   template< int R >
   struct real
	 : pegtl::sor< pegtl::seq< sign, ureal< R > >,
		       pegtl::seq< pegtl::one< '+' >, naninf >,
		       pegtl::seq< pegtl::one< '-' >, naninf > > {};

   template< int R >
   struct complex
	 : pegtl::sor< real< R >,
		       pegtl::seq< real< R >, pegtl::one< '@' >, real< R > >,
		       pegtl::seq< real< R >, pegtl::one< '+' >, ureal< R >, pegtl::one< i > >,
		       pegtl::seq< real< R >, pegtl::one< '-' >, ureal< R >, pegtl::one< i > >,
		       pegtl::seq< real< R >, pegtl::one< '+' >, naninf, pegtl::one< i > >,
		       pegtl::seq< real< R >, pegtl::one< '-' >, naninf, pegtl::one< i > >,
		       pegtl::seq< real< R >, pegtl::one< '+' >, pegtl::one< i > >,
		       pegtl::seq< real< R >, pegtl::one< '-' >, pegtl::one< i > >,
		       pegtl::seq< pegtl::one< '+' >, ureal< R >, pegtl::one< i > >,
		       pegtl::seq< pegtl::one< '-' >, ureal< R >, pegtl::one< i > >,
		       pegtl::seq< pegtl::one< '+' >, naninf, pegtl::one< i > >,
		       pegtl::seq< pegtl::one< '-' >, naninf, pegtl::one< i > >,
		       pegtl::seq< pegtl::one< '+' >, pegtl::one< 'i' > >,
		       pegtl::seq< pegtl::one< '-' >, pegtl::one< 'i' > > > {};

   template< int R >
   struct num
	 : pegtl::seq< prefix< R >, complex< R > > {};

   typedef pegtl::sor< num< 2 >,
		       num< 8 >,
		       num< 16 >,
		       num< 10 > > number_;

   struct number
	 : padded_lexeme< number_ > {};

   struct line_ending
	 : pegtl::eol {};

   struct hex_scalar_value
	 : pegtl::plus< hex_digit > {};

   struct character_tabulation
	 : pegtl::blank {};

   struct intraline_whitespace
	 : character_tabulation {};

   struct special_initial
	 : pegtl::one< '!', '$', '%', '&', '*', '/', ':', '<', '=', '>', '?', '^', '_', '~' > {};

   struct letter
	 : pegtl::alpha {};

   struct constituent
	 : letter {};

   struct whitespace
	 : pegtl::space {};

   struct delimiter
	 : pegtl::sor< pegtl::one< '(', ')', '[', ']', '"', ';', '#' >, whitespace, pegtl::eol > {};

   struct inline_hex_escape
	 : pegtl::seq< pegtl::string< '\\', x >, pegtl::seq< hex_scalar_value, pegtl::one< ';' > > > {};

   struct string_element
	 : pegtl::sor< pegtl::not_one< '"', '\\' >, pegtl::seq< pegtl::at_one< '\\' >, pegtl::sor< pegtl::seq< pegtl::one< '\\' >, pegtl::one< a, b, t, n, v, f, r, '"', '\\' > >, pegtl::seq< pegtl::one< '\\' >, intraline_whitespace, line_ending, intraline_whitespace >, inline_hex_escape > > > {};

   typedef pegtl::seq< pegtl::one< '"' >, pegtl::seq< pegtl::star< string_element >, pegtl::one< '"' > > > string_;

   struct string
	 : padded_lexeme< string_ > {};

   struct character_name
	 : pegtl::sor< pegtl::string< n, u, l >, pegtl::string< a, l, a, r, m >, pegtl::string< b, a, c, k, s, p, a, c, e >, pegtl::string< t, a, b >, pegtl::string< l, i, n, e, f, e, e, d >, pegtl::string< n, e, w, l, i, n, e >, pegtl::string< v, t, a, b >, pegtl::string< p, a, g, e >, pegtl::string< r, e, t, u, r, n >, pegtl::string< e, s, c >, pegtl::string< s, p, a, c, e >, pegtl::string< d, e, l, e, t, e > > {};

   typedef pegtl::seq< pegtl::seq< pegtl::string< '#', '\\' >, pegtl::sor< pegtl::seq< pegtl::one< x >, hex_scalar_value >, character_name, pegtl::any > >, pegtl::at< delimiter > > character_;

   struct character
	 : padded_lexeme< character_ > {};

   typedef pegtl::seq< pegtl::seq< pegtl::one< '#' >, pegtl::one< t, 'T', f, 'F' > >, pegtl::at< delimiter > > boolean_;

   struct boolean
	 : padded_lexeme< boolean_ > {};

   struct subsequent;

   struct peculiar_identifier
	 : pegtl::seq< pegtl::sor< pegtl::one< '+' >, pegtl::one< '-' >, pegtl::string< '.', '.', '.' >, pegtl::seq< pegtl::string< '-', '>' >, pegtl::star< subsequent > > >, pegtl::at< delimiter > > {};

   struct special_subsequent
	 : pegtl::one< '+', '-', '.', '@' > {};

   struct initial
	 : pegtl::sor< constituent, special_initial, inline_hex_escape > {};

   struct subsequent
	 : pegtl::sor< initial, digit< 10 >, special_subsequent > {};

   typedef pegtl::seq< pegtl::sor< pegtl::seq< initial, pegtl::star< subsequent > >, peculiar_identifier >, pegtl::at< delimiter > > identifier_;

   struct identifier
	 : padded_lexeme< identifier_ > {};

   struct whitespace;
   struct comment;

   // We don't need a star here because interlexeme_space is used only as second argument to pegtl::pad, which automatically adds a star.

   struct interlexeme_space
	 : pegtl::sor< whitespace, comment > {};

   struct comment_text;
   struct nested_comment;

   struct comment_cont
	 : pegtl::seq< nested_comment, comment_text > {};

   struct comment_text
	 : pegtl::until< pegtl::sor< pegtl::string< '#', '|' >, pegtl::string< '|', '#' > > > {};

   struct nested_comment
	 : pegtl::seq< pegtl::string< '#', '|' >, pegtl::seq< comment_text, pegtl::star< comment_cont >, pegtl::string< '|', '#' > > > {};

   struct datum;

   struct comment
	 : pegtl::sor< pegtl::seq< pegtl::one< ';' >, pegtl::until< pegtl::eol > >, nested_comment, pegtl::seq< pegtl::string< '#', ';' >, interlexeme_space, datum >, pegtl::string< '#', '!', r, '6', r, s > > {};

   struct u8
	 : number {};

   struct bytevector
	 : pegtl::seq< padded_lexeme< pegtl::string< '#', v, u, '8', '(' > >, pegtl::seq< pegtl::star< u8 >, padded_lexeme< pegtl::one< ')' > > > > {};

   struct vector
	 : pegtl::seq< padded_lexeme< pegtl::string< '#', '(' > >, pegtl::seq< pegtl::star< datum >, padded_lexeme< pegtl::one< ')' > > > > {};

   struct abbrev_prefix
	 : padded_lexeme< pegtl::sor< pegtl::one< '`' >, pegtl::one< '\'' >, pegtl::one< ',' >, pegtl::string< ',', '@' >, pegtl::string< '#', '`' >, pegtl::string< '#', '\'' >, pegtl::string< '#', ',' >, pegtl::string< '#', ',', '@' > > > {};

   struct abbreviation
	 : pegtl::seq< abbrev_prefix, datum > {};

   struct datum;

   struct list
	 : pegtl::sor< pegtl::seq< padded_lexeme< pegtl::one< '(' > >, pegtl::seq< pegtl::sor< pegtl::seq< pegtl::plus< datum >, pegtl::ifmust< pegtl::one< '.' >, pegtl::at< delimiter > >, datum >, pegtl::star< datum > >, padded_lexeme< pegtl::one< ')' > > > >,
		       pegtl::seq< padded_lexeme< pegtl::one< '[' > >, pegtl::seq< pegtl::sor< pegtl::seq< pegtl::plus< datum >, pegtl::ifmust< pegtl::one< '.' >, pegtl::at< delimiter > >, datum >, pegtl::star< datum > >, padded_lexeme< pegtl::one< ']' > > > >,
		       abbreviation > {};

   struct compound_datum
	 : pegtl::sor< bytevector, vector, list > {};

   struct symbol
	 : identifier {};

   struct lexeme_datum
	 : pegtl::sor< boolean, number, character, string, symbol > {};

   struct datum
	 : pegtl::sor< lexeme_datum, compound_datum > {};

   struct r6rs
	 : pegtl::until< pegtl::space_until_eof, datum > {};

} // scheme

int main( int argc, char ** argv )
{
   if ( argc == 1 ) {
      pegtl::print_rules< scheme::datum >();
   }
   for ( int arg = 1; arg < argc; ++arg ) {
      pegtl::smart_parse_file< scheme::r6rs >( false, argv[ arg ] );
      PEGTL_PRINT( "input from file " << argv[ arg ] << " accepted" );
   }
   return 0;
}
