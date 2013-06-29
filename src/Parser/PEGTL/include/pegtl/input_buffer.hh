// Copyright (c) 2008 Dr. Colin Hirsch 
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_INPUT_BUFFER_HH
#define COHI_PEGTL_INPUT_BUFFER_HH


namespace pegtl
{
   // Classes for input given as range of input iterators
   // with automatic minimal buffering for back-tracking.

   template< typename Iterator >
   class buffer_impl : private nocopy< buffer_impl< Iterator > >
   {
   public:
      buffer_impl( const Iterator begin, const Iterator end )
	    : m_record( false ),
	      m_count( 0 ),
	      m_runpos( 0 ),
	      m_vecpos( 0 ),
	      m_run( begin ),
	      m_end( end )
      { }

      ~buffer_impl()
      {
	 PEGTL_PRINT( "pegtl: buffer runpos " << m_runpos << " vecpos " << m_vecpos << " capacity " << m_vector.capacity() );
      }

      typedef typename std::iterator_traits< Iterator >::value_type value_type;

      bool eof( const size_t offset )
      {
	 //	 PEGTL_DEBUG( offset );
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );
	 if ( offset > m_runpos ) {
	    if ( m_record ) {
	       record_true( offset );
	    }
	    else {
	       record_false( offset );
	    }
	 }
	 assert( offset <= m_runpos );
	 return ( offset == m_runpos ) && ( m_run == m_end );
      }

      bool record() const
      {
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );
	 return m_record;
      }

      void record_start( const size_t offset )
      {
	 //	 ++m_count;
	 //	 PEGTL_DEBUG( offset );
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );

	 if ( m_vecpos != offset ) {
	    m_vector.clear();
	    m_vecpos = offset;
	 }
	 m_record = true;
      }

      void record_stop()
      {
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );
	 m_record = false;
      }

      void record_true( const size_t offset )
      {
	 while ( offset > m_runpos ) {
	    assert( m_run != m_end );
	    m_vector.push_back( * m_run );
	    ++m_run;
	    ++m_runpos;
	 }
      }

      void record_false( const size_t offset )
      {
	 m_vector.clear();
	 m_vecpos = m_runpos;
	 record_true( offset );
      }

      const value_type & operator[] ( const size_t offset )
      {
	 //	 PEGTL_DEBUG( offset );
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );

	 if ( m_record ) {
	    record_true( offset );
	 }
	 while ( offset > m_runpos ) {
	    assert( m_run != m_end );
	    ++m_run;
	    ++m_runpos;
	    //	    PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() );
	 }
	 if ( ( offset >= m_vecpos ) && ( offset < m_vecpos + m_vector.size() ) ) {
	    //	    PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() );
	    //	    PEGTL_DEBUG( m_runpos << " C " << int( m_vector[ offset - m_vecpos ] ) << " " << m_vector[ offset - m_vecpos ] );
	    return m_vector[ offset - m_vecpos ];
	 }
	 assert( ! eof( offset ) );
	 if ( offset == m_runpos ) {
	    //	    PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() );
	    //	    PEGTL_DEBUG( m_runpos << " B " << int( * m_run ) << " '" << ( * m_run ) << "'" );
	    return * m_run;
	 }
	 //	 PEGTL_DEBUG( m_record << ( m_run == m_end ) << " " << m_runpos << " " << m_vecpos << " " << m_vector.size() << " " << m_count );
	 PEGTL_THROW( "buffer iterator access error" );
      }

   private:
      bool m_record;
      size_t m_count;

      size_t m_runpos;
      size_t m_vecpos;

      Iterator m_run;
      const Iterator m_end;

      std::vector< value_type > m_vector;
   };

   template< bool Must, typename Input > class buffer_marker;

   template< typename Input >
   class buffer_marker< true, Input >
   {
   public:
      explicit
      buffer_marker( Input & )
      { }

      bool operator() ( const bool success )
      {
	 return success;
      }
   };

   template< typename Input >
   class buffer_marker< false, Input >
   {
   public:
      explicit
      buffer_marker( Input & in )
	    : m_input( & in ),
	      m_stop( in.record_start() ),
	      m_iterator( in.here() )
      { }

      ~buffer_marker()
      {
	 if ( m_input ) {
	    m_input->record_stop( m_stop );
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
	    m_input->record_stop( m_stop );
	    m_input = 0;
	 }
	 return success;
      }

   private:
      Input * m_input;
      const bool m_stop;
      const typename Input::iterator m_iterator;
   };

   template< typename Iterator, typename Location >
   class buffer_iterator : public std::iterator< std::forward_iterator_tag, typename std::iterator_traits< Iterator >::value_type >
   {
   public:
      buffer_iterator( const size_t offset, buffer_impl< Iterator > & buffer )
	    : m_offset( offset ),
	      m_buffer( & buffer )
      { }

      Location location() const
      {
	 return Location( m_offset );
      }

      typedef typename std::iterator_traits< Iterator >::value_type value_type;

      value_type operator* () const
      {
	 return ( * m_buffer )[ m_offset ];
      }

      const value_type * operator-> () const
      {
	 return & operator* ();
      }

      buffer_iterator operator++ ( int )
      {
	 buffer_iterator nrv( * this );
	 this->operator++ ();
	 return nrv;
      }

      const buffer_iterator & operator++ ()
      {
	 ++m_offset;
	 return * this;
      }

      bool operator== ( const buffer_iterator & i ) const
      {
	 return m_offset == i.m_offset;
      }

      bool operator!= ( const buffer_iterator & i ) const
      {
	 return m_offset == i.m_offset;
      }

   public:
      size_t internal_offset() const
      {
	 return m_offset;
      }

   private:
      size_t m_offset;
      buffer_impl< Iterator > * m_buffer;
   };

   template< typename Iterator, typename Location >
   class buffer_input : private nocopy< buffer_input< Iterator, Location > >
   {
   public:
      buffer_input( const Iterator begin, const Iterator end )
	    : m_buffer( begin, end ),
	      m_run( 0, m_buffer )
      { }

      typedef Location location_type;
      typedef buffer_iterator< Iterator, Location > iterator;
      typedef typename std::iterator_traits< Iterator >::value_type value_type;

      template< bool Must >
      struct marker : public buffer_marker< Must, buffer_input >
      {
	 explicit
	 marker( buffer_input & input )
	       : buffer_marker< Must, buffer_input >( input )
	 { }
      };

      bool eof() const
      {
	 return m_buffer.eof( m_run.internal_offset() );
      }

      bool record_start()
      {
	 if ( m_buffer.record() ) {
	    return false;
	 }
	 m_buffer.record_start( m_run.internal_offset() );
	 return true;
      }

      void record_stop( const bool stop )
      {
	 if ( stop ) {
	    m_buffer.record_stop();
	 }
      }

      const iterator & here() const
      {
	 return m_run;
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

   protected:
      mutable buffer_impl< Iterator > m_buffer;
      iterator m_run;

      void throw_at_eof() const
      {
	 if ( eof() ) {
	    PEGTL_THROW( "attempt to read beyond end of input" );
	 }
      }
   };

} // pegtl

#endif
