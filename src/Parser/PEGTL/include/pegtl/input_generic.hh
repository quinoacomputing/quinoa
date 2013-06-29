// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_INPUT_GENERIC_HH
#define COHI_PEGTL_INPUT_GENERIC_HH


namespace pegtl
{
   // The following classes are helper classes to simplify the implementation
   // of rules: used correctly, they make adherence to 'do not consume input
   // on failure' easy.
   // They are modelled after a transactional style of programming.
   // The constructor is the 'begin', and 'operator()' when called with
   // a value of 'true' is the 'commit'. A 'rollback' is performed either
   // when 'operator()' is called with 'false', or when the destructor is
   // called without 'operator()' having been called earlier; this provides
   // for correct behaviour in the presence of exceptions (side note: in C++,
   // the destructor is the most important 'thing' for exception safety, not
   // the try-catch block, as in some other languages I do not want to name).

   // Class character acts as proxy to the current character in the input;
   // a 'commit' consumes the character (thereby moving the input to the next
   // position), a 'rollback' is a nop.

   template< typename Input >
   struct character
   {
   public:
      explicit
      character( Input & in )
	    : m_input( in ),
	      m_value( in.peek() )
      { }

      typedef typename Input::value_type value_type;

      operator const value_type & () const
      {
	 return m_value;
      }

      bool operator() ( const bool success ) const
      {
	 if ( success ) {
	    m_input.bump();
	 }
	 return success;
      }

   private:
      Input & m_input;
      const value_type m_value;
   };

   // Please see the supplied rule classes for usage examples...

   class dummy_location
   {
   public:
      int operator() ( const int c )
      {
	 return c;
      }

      void write_to( std::ostream & o ) const
      {
	 o << '?';
      }
   };

   class offset_location
   {
   public:
      explicit
      offset_location( const size_t offset = 0 )
	    : m_offset( offset )
      { }

      int operator() ( const int c )
      {
	 ++m_offset;
	 return c;
      }

      void write_to( std::ostream & o ) const
      {
	 o << m_offset;
      }

   private:
      size_t m_offset;
   };

   inline std::ostream & operator<< ( std::ostream & o, const dummy_location & w )
   {
      w.write_to( o );
      return o;
   }

   inline std::ostream & operator<< ( std::ostream & o, const offset_location & w )
   {
      w.write_to( o );
      return o;
   }

} // pegtl

#endif
