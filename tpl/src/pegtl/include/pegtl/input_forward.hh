// Copyright (c) 2008 Dr. Colin Hirsch 
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_INPUT_FORWARD_HH
#define COHI_PEGTL_INPUT_FORWARD_HH


namespace pegtl
{
   // Classes for input given as range of forward iterators.

   // Class marker remembers (marks) the current position in the input; a
   // 'commit' is a nop, a 'rollback' rewinds to the initial position; the
   // 'Must' flag is used to optimise: if a rule must succeed (or else an
   // exception aborts the parser), then we don't need to enable this kind
   // of backtracking.

   template< bool Must, typename Input > class forward_marker;

   template< typename Input >
   class forward_marker< true, Input >
   {
   public:
      explicit
      forward_marker( Input & )
      { }

      bool operator() ( const bool success )
      {
	 return success;
      }
   };

   template< typename Input >
   class forward_marker< false, Input >
   {
   public:
      explicit
      forward_marker( Input & in )
	    : m_input( & in ),
	      m_iterator( in.here() )
      { }

      ~forward_marker()
      {
	 if ( m_input ) {
	    m_input->jump( m_iterator );
	 }
      }

      typedef typename Input::iterator iterator;

      iterator here() const
      {
	 return m_iterator;
      }

      bool operator() ( const bool success )
      {
	 if ( success ) {
	    m_input = 0;
	 }
	 return success;
      }

   private:
      Input * m_input;
      const typename Input::iterator m_iterator;
   };

   template< typename Iterator, typename Location >
   class forward_iterator : public std::iterator< std::forward_iterator_tag, typename std::iterator_traits< Iterator >::value_type >
   {
   public:
      forward_iterator()
      { }

      explicit
      forward_iterator( const Iterator & it )
	    : m_iterator( it )
      { }

      typedef typename std::iterator_traits< Iterator >::value_type value_type;

      value_type operator* () const
      {
	 return * m_iterator;
      }

      const value_type * operator-> () const
      {
	 return m_iterator.operator->();
      }

      forward_iterator operator++ ( int )
      {
	 forward_iterator nrv( * this );
	 this->operator++ ();
	 return nrv;
      }

      const forward_iterator & operator++ ()
      {
	 m_location( * m_iterator );
	 ++m_iterator;
	 return * this;
      }

      const Location & location() const
      {
	 return m_location;
      }

      bool operator== ( const forward_iterator & i ) const
      {
	 return m_iterator == i.m_iterator;
      }

      bool operator!= ( const forward_iterator & i ) const
      {
	 return m_iterator != i.m_iterator;
      }

      void assign( const Iterator & it )
      {
	 m_iterator = it;
	 m_location = Location();
      }

   private:
      Location m_location;
      Iterator m_iterator;
   };

   template< typename Iterator, typename Location >
   class forward_input : private nocopy< forward_input< Iterator, Location > >
   {
   public:
      forward_input( const Iterator begin, const Iterator end )
	    : m_run( begin ),
	      m_begin( begin ),
	      m_end( end )
      { }

      typedef Location location_type;
      typedef forward_iterator< Iterator, Location > iterator;
      typedef typename std::iterator_traits< Iterator >::value_type value_type;

      template< bool Must >
      struct marker : public forward_marker< Must, forward_input >
      {
	 explicit
	 marker( forward_input & input )
	       : forward_marker< Must, forward_input >( input )
	 { }
      };

      bool eof() const
      {
	 return m_run == m_end;
      }

      const iterator & here() const
      {
	 return m_run;
      }

      forward_input & rewind()
      {
	 m_run = m_begin;
	 return * this;
      }

      Location location() const
      {
	 return m_run.location();
      }

      void bump()
      {
	 throw_at_eof();
	 ++m_run;
      }

      value_type peek() const
      {
	 throw_at_eof();
	 return * m_run;
      }

      void jump( const iterator iter )
      {
	 m_run = iter;
      }

      std::string debug_escape( const iterator & begin, const iterator & end ) const
      {
	 std::string nrv;

	 for ( iterator run = begin; run != end; ++run ) {
	    escape_impl( nrv, * run );
	 }
	 return nrv;
      }

   protected:
      iterator m_run;

      const iterator m_begin;
      const iterator m_end;

      void throw_at_eof() const
      {
	 if ( eof() ) {
	    PEGTL_THROW( "attempt to read beyond end of input" );
	 }
      }
   };

} // pegtl

#endif
