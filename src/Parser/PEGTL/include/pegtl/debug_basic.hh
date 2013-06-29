// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_DEBUG_BASIC_HH
#define COHI_PEGTL_DEBUG_BASIC_HH


namespace pegtl
{
   template< typename Rule, typename Input, typename Debug >
   struct basic_guard
	 : private nocopy< basic_guard< Rule, Input, Debug > >
   {
      typedef typename Input::location_type Location;

      basic_guard( Location && w, counter & t )
	    : m_location( std::move( w ) ),
	      m_counter( t )
      {
	 m_counter.enter();
      }

      ~basic_guard()
      {
	 m_counter.leave();
      }

      bool operator() ( const bool result, const bool must ) const
      {
	 if ( ( ! result ) && must ) {
	    PEGTL_THROW( "parsing aborted at " << m_location );
	 }
	 return result;
      }

      const Location & location() const
      {
	 return m_location;
      }

   protected:
      const Location m_location;
      counter & m_counter;
   };


   struct basic_debug : public debug_base
   {
      template< typename TopRule >
      explicit
      basic_debug( const tag< TopRule > & help )
	    : m_printer( help )
      { }

      template< bool Must, typename Rule, typename Input, typename ... States >
      bool match( Input & in, States && ... st )
      {
	 const basic_guard< Rule, Input, basic_debug > d( in.location(), m_counter );

	 try {
	    return d( Rule::template match< Must >( in, * this, std::forward< States >( st ) ... ), Must );
	 }
	 catch ( const parse_error & ) {
	    PEGTL_LOGGER( * this, "#" << m_counter.nest() << " @" << d.location() << ": " << m_printer.template rule< Rule >() );
	    throw;
	 }
      }

   protected:
      counter m_counter;
      printer m_printer;
   };

} // pegtl

#endif
