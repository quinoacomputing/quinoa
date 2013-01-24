// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_PARSE_FILENAME_HH
#define COHI_PEGTL_PARSE_FILENAME_HH

namespace pegtl
{
   template< typename Location >
   struct file_input
	 : private file_mapper,
	   public forward_input< file_mapper::iterator, Location >
   {
      explicit
      file_input( const std::string & filename )
	    : file_mapper( filename ),
	      forward_input< file_mapper::iterator, Location >( file_mapper::begin(), file_mapper::end() )
      { }

      typedef typename forward_input< file_mapper::iterator, Location >::iterator iterator;
   };

   typedef file_input< ascii_location > ascii_file_input;
   typedef file_input< dummy_location > dummy_file_input;

   // Functions to parse input given the filename as std::string.

   // Wrapper functions that add another convenience layer: instantiation
   // and initialisation of the input class. See file parse_generic.hh for
   // the wrapped functions.

   // The *_parse_file_* functions set up the parser to parse a file; the
   // functions here use mmap(2), alternatively the function read_string()
   // defined in utilities.hh can be used to read a file into a std::string
   // (which can then be parsed by one of the parse_string functions).

   template< typename TopRule, typename ... States >
   void dummy_parse_file( const std::string & filename, States && ... st )
   {
      dummy_file_input in( filename );
      dummy_parse< TopRule >( in, std::forward< States >( st ) ... );
   }

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void basic_parse_file( const std::string & filename, States && ... st )
   {
      file_input< Location > in( filename );
      basic_parse< TopRule >( in, std::forward< States >( st ) ... );
   }

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void trace_parse_file( const bool trace, const std::string & filename, States && ... st )
   {
      file_input< Location > in( filename );
      trace_parse< TopRule >( trace, in, std::forward< States >( st ) ... );
   }

   // Please read the comment on the smart_parse_* functions in parse_generic.hh!

   template< typename TopRule, typename Location = ascii_location, typename ... States >
   void smart_parse_file( const bool trace, const std::string & filename, States && ... st )
   {
      file_input< Location > in( filename );
      smart_parse< TopRule >( trace, in, std::forward< States >( st ) ... );
   }

} // pegtl

#endif
