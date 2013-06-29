// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

// Include the only public header file for the PEGTL.

#include <pegtl.hh>

namespace example
{
   using namespace pegtl;

   // Define a new grammar, i.e. a new rule, combining
   // some of the rules that are part of the PEGTL.

   struct grammar
	 : seq< alpha, until< eol, sor< alpha, digit > > > {};

   // The atomic rule 'alpha' matches any upper- or lower-case ASCII character.
   // The atomic rule 'digit' matches any ASCII digit.
   // The atomic rule 'eol' matches any ASCII end-of-line, and it ALSO matches end-of-file.

   // The rule combinator 'seq' stands for concatenation
   // The rule combinator 'sor' stands for ordered-choice, or sequenced-or.
   // The rule combinator 'until' matches its second argument rule until the first argument rule matches.

   // Put together, this grammar is equivalent to the regular expression '^[[:alpha:]]([[:alpha:]]|[[:digit:]])*$'.

   // The initial '^' is implicit; matching a rule never "just skips" input when the rule does not match.
   // The final '$' is explicit in form of the eol-rule; otherwise matching could successfully finish earlier.

   // See the PEGTL documentation on why grammar is defined as a struct, rather than a simple typedef.

} // example

int main( int argc, char ** argv )
{
   // Parse all command-line arguments with the grammar.

   for ( int i = 1; i < argc; ++i ) {
      // The member of the parse-family of functions used here
      // - uses a basic_debug for diagnostics,
      // - takes its input from a std::string,
      // and is therefore named basic_parse_string().
      // As all parse functions, it returns on success, and
      // throws a pegtl::parse_error on failure.

      pegtl::basic_parse_string< example::grammar >( argv[ i ] );
   }
   return 0;
}
