// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_DEBUG_TRACE_HH
#define COHI_PEGTL_DEBUG_TRACE_HH


namespace pegtl
{
   template< typename Rule, typename Input, typename Debug >
   struct trace_guard
	 : public basic_guard< Rule, Input, Debug >
   {
   public:
      trace_guard( Input & in, Debug & de, counter & t, printer & p )
	    : basic_guard< Rule, Input, Debug >( in.location(), t ),
	      m_input( in ),
	      m_debug( de ),
	      m_printer( p ),
	      m_begin( in.here() ),
	      m_result( 2 ),
	      m_rule_counter( t.rule() )
      {
	 trace( "start  " );
      }

      ~trace_guard()
      {
	 const char * const msgs[] = { "failure", "success", "unwind " };
	 trace( msgs[ m_result ] );
      }

      bool operator() ( const bool result, const bool must )
      {
	 return basic_guard< Rule, Input, Debug >::operator() ( ( m_result = result ), must );
      }

   private:
      Input & m_input;
      Debug & m_debug;
      printer & m_printer;
      const typename Input::iterator m_begin;
      unsigned m_result;
      const unsigned m_rule_counter;

      void trace( const char * msg )
      {
	 PEGTL_LOGGER( m_debug, msg << " flags " << m_result << " rule " << std::setw( 4 ) << m_rule_counter << " nest " << std::setw( 3 ) << this->m_counter.nest() << " at " << m_input.location() << " expression " << this->m_printer.template rule< Rule >() << " input \"" << m_input.debug_escape( m_begin, m_input.here() ) << "\"" );
      }
   };

   struct trace_debug : public debug_base
   {
      template< typename TopRule >
      explicit
      trace_debug( const tag< TopRule > & help, const bool trace = true )
	    : m_trace( trace ),
	      m_printer( help )
      { }

      bool trace() const
      {
	 return m_trace;
      }

      void set_trace( const bool trace )
      {
	 m_trace = trace;
      }

      template< bool Must, typename Rule, typename Input, typename ... States >
      bool match( Input & in, States && ... st )
      {
	 return m_trace ? match_trace< Must, Rule >( in, std::forward< States >( st ) ... ) : match_basic< Must, Rule >( in, std::forward< States >( st ) ... );
      }

   protected:
      bool m_trace;
      counter m_counter;
      printer m_printer;

   private:
      template< bool Must, typename Rule, typename Input, typename ... States >
      bool match_basic( Input & in, States && ... st )
      {
	 basic_guard< Rule, Input, trace_debug > d( in.location(), m_counter );

	 try {
	    return d( Rule::template match< Must >( in, * this, std::forward< States >( st ) ... ), Must );
	 }
	 catch ( const parse_error & ) {
	    PEGTL_LOGGER( * this, "#" << m_counter.nest() << " @" << d.location() << ": " << m_printer.template rule< Rule >() );
	    throw;
	 }
      }

      template< bool Must, typename Rule, typename Input, typename ... States >
      bool match_trace( Input & in, States && ... st )
      {
	 trace_guard< Rule, Input, trace_debug > d( in, * this, m_counter, m_printer );

	 try {
	    return d( Rule::template match< Must >( in, * this, std::forward< States >( st ) ... ), Must );
	 }
	 catch ( const parse_error & ) {
	    PEGTL_LOGGER( * this, "#" << m_counter.nest() << " @" << d.location() << ": " << m_printer.template rule< Rule >() );
	    throw;
	 }
      }

   };

} // pegtl

#endif
