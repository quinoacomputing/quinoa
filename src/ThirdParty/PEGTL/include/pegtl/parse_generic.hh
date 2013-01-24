// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_PARSE_GENERIC_HH
#define COHI_PEGTL_PARSE_GENERIC_HH

namespace pegtl
{
   // The following function is the central wrapper used for parsing.

   template< typename TopRule, typename Input, typename Debug, typename ... States >
   void parse( Input & in, Debug & de, States && ... st )
   {
      // Unless Debug is broken, the return value is always true [because the Debug
      // must convert a local failure (return false) into a global failure (exception)
      // when Must == true, as is the case here].
      const bool b = de.template match< true, TopRule >( in, std::forward< States >( st ) ... );
      assert( b );
   }

   // The next functions add one convenience layer: they take care of instantiating
   // one (or more) of the debug classes supplied with the library, corresponding to
   // their respective names.

   // Parse with a dummy_debug (full speed, no diagnostics).

   template< typename TopRule, typename Input, typename ... States >
   void dummy_parse( Input & in, States && ... st )
   {
      dummy_debug de;
      parse< TopRule >( in, de, std::forward< States >( st ) ... );
   }

   // Parse with a dummy_debug (slower, diagnostics on error).

   template< typename TopRule, typename Input, typename ... States >
   void basic_parse( Input & in, States && ... st )
   {
      basic_debug de( tag< TopRule >( 0 ) );
      parse< TopRule >( in, de, std::forward< States >( st ) ... );
   }

   // Parse with a trace_debug (fixed overhead depending on the grammar,
   // diagnostics on error, plus optionally trace of all rule invocations).
   // The first bool argument called trace enables or disables rule tracing.

   template< typename TopRule, typename Input, typename ... States >
   void trace_parse( const bool trace, Input & in, States && ... st )
   {
      trace_debug de( tag< TopRule >( 0 ), trace );
      parse< TopRule >( in, de, std::forward< States >( st ) ... );
   }

   // Parse with a dummy_debug and, when that fails, use a trace_debug, i.e.
   // try the fastest way, with neither diagnostics nor safety-nets first,
   // and fall-back to a more verbose version only when necessary. In case
   // of errors, the whole input is parsed two times. Note that the state
   // is not cleared in between the two parsing runs; it might be necessary
   // for the application to copy'n'paste this function with an appropriate
   // preparation of the state before starting the second run...!

   template< typename TopRule, typename Input, typename ... States >
   void smart_parse( const bool trace, Input & in, States && ... st )
   {
      try {
	 dummy_parse< TopRule >( in, std::forward< States >( st ) ... );
      }
      catch ( const parse_error & ) {
	 trace_parse< TopRule >( trace, in.rewind(), std::forward< States >( st ) ... );
      }
   }

} // pegtl

#endif
