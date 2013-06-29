// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace example
{
   using namespace pegtl;

   // Class capture is both a rule and an action.

   // As action, it takes the matched sub-string passed to its apply() method
   // and stores it in the first state argument with the template argument as
   // key, i.e. for state m, matched sub-string s, and template argument Key
   // the assignment m[ Key ] = s is performed.

   // As rule, it retrieves m[ Key ] and matches the input against the resulting
   // string; it effectively behaves like the string<> rule, but with the string
   // to match against fixed at run-time by the map in the state. When no entry
   // can be found in the map for the given key the rule fails.

   // In order to make sense, the grammar must, for every given key, ensure that
   // capture is first used as action to store a string in the map, and only later
   // as rule that retrieves and uses the previously stored string.

   // This example matches all strings that consist of the same two sequences of
   // digits, separated by a non-empty sequence of tabs and spaces.

   struct grammar
	 : seq< ifapply< plus< digit >, capture< 42 > >, plus< blank >, capture< 42 >, eof > {};

} // example

int main( int argc, char ** argv )
{
   pegtl::capture_map map;

   for ( int i = 1; i < argc; ++i ) {
      pegtl::basic_parse_string< example::grammar >( argv[ i ], map );
      PEGTL_PRINT( "input " << argv[ i ] << " map entry " << map[ 42 ] );
   }
   return 0;
}
