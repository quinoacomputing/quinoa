// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace grammar
{
   // The first non-trivial grammar used during development and debugging of the
   // library. This grammar recognises a small subset of parsing expressions.

   using namespace pegtl;

   struct read_expr;

   struct read_comment
	 : seq< one< '#' >, until< eol > > {};

   struct read_terminal_char
	 : ifmust< one< '\'' >, seq< not_one< '\'' >, one< '\'' > > > {};

   struct read_terminal_string
	 : ifmust< one< '"' >, seq< star< not_one< '"' > >, one< '"' > > > {};

   struct read_infix
	 : pad< one< '/', '.' >, blank > {};

   struct read_prefix
	 : pad< one< '!', '&' >, blank > {};

   struct read_postfix
	 : pad< one< '+', '*', '?' >, blank > {};

   struct read_paren
	 : ifmust< pad_one< '(', blank >, seq< read_expr, pad_one< ')', blank > > > {};

   struct read_atomic
	 : sor< read_terminal_char, read_terminal_string > {};

   struct read_body
	 : sor< read_paren, read_atomic > {};

   struct read_expr
	 : seq< opt< read_prefix >, read_body, opt< read_postfix >, ifthen< read_infix, read_expr > > {};

   struct read_rule
	 : seq< identifier, pad_one< '=', blank >, read_expr > {};

   struct read_line
	 : sor< read_comment, read_rule > {};

   struct read_file
	 : until< space_until_eof, read_line > {};

} // grammar

int main( int argc, char ** argv )
{
   for ( int arg = 1; arg < argc; ++arg ) {
      pegtl::basic_parse_string< grammar::read_file >( argv[ arg ] );
      std::cerr << "input " << argv[ arg ] << " valid\n";
   }
   return 0;
}
