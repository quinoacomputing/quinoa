// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_UTILITIES_HH
#define COHI_PEGTL_UTILITIES_HH

namespace pegtl
{
   template< typename T >
   struct tag
   {
      typedef T type;

      tag()
      { }

      explicit
      tag( const int )
      { }

      explicit
      tag( const volatile T * )
      { }
   };

#define PEGTL_PRINT( MeSSaGe )				\
   do { std::cerr << MeSSaGe << '\n'; } while( 0 )

#define PEGTL_TRACE							\
   PEGTL_PRINT( __FILE__ << " : " << __LINE__ << " : " << __PRETTY_FUNCTION__ )

#define PEGTL_DEBUG( MeSSaGe )						\
   PEGTL_PRINT( __FILE__ << " : " << __LINE__ << " : " << MeSSaGe )

   struct parse_error
	 : public std::runtime_error
   {
      explicit
      parse_error( const std::string & message )
	    : std::runtime_error( message )
      { }
   };

#define PEGTL_THROW( MeSSaGe )						\
   do { std::ostringstream oss; oss << MeSSaGe; throw ::pegtl::parse_error( oss.str() ); } while( 1 )

#define PEGTL_LOGGER( DeBuG, MeSSaGe )			\
   do { std::ostringstream oss; oss << MeSSaGe; ( DeBuG ).log( oss.str() ); } while( 0 )

   inline void escape_impl( std::string & result, const int i )
   {
      switch ( i )
      {
	 case '"':
	    result += "\\\"";
	    break;
	 case '\\':
	    result += "\\\\";
	    break;
	 case '\a':
	    result += "\\a";
	    break;
	 case '\b':
	    result += "\\b";
	    break;
	 case '\t':
	    result += "\\t";
	    break;
	 case '\n':
	    result += "\\n";
	    break;
	 case '\r':
	    result += "\\r";
	    break;
	 case '\v':
	    result += "\\v";
	    break;
	 case 32: case 33:          case 35: case 36: case 37: case 38: case 39:
         case 40: case 41: case 42: case 43: case 44: case 45: case 46: case 47: case 48: case 49:
         case 50: case 51: case 52: case 53: case 54: case 55: case 56: case 57: case 58: case 59:
         case 60: case 61: case 62: case 63: case 64: case 65: case 66: case 67: case 68: case 69:
         case 70: case 71: case 72: case 73: case 74: case 75: case 76: case 77: case 78: case 79:
         case 80: case 81: case 82: case 83: case 84: case 85: case 86: case 87: case 88: case 89:
         case 90: case 91:          case 93: case 94: case 95: case 96: case 97: case 98: case 99:
         case 100: case 101: case 102: case 103: case 104: case 105: case 106: case 107: case 108: case 109:
         case 110: case 111: case 112: case 113: case 114: case 115: case 116: case 117: case 118: case 119:
         case 120: case 121: case 122: case 123: case 124: case 125: case 126:
	    result += char( i );
	    break;
	 default: {
	    char tmp[ 12 ];
	    ::snprintf( tmp, sizeof( tmp ), "\\u%04x", i );
	    result += tmp;
	 }  break;
      }
   }

   inline std::string escape( const int i )
   {
      std::string nrv;
      escape_impl( nrv, i );
      return nrv;
   }

   inline std::string escape( const std::string & s )
   {
      std::string nrv;

      for ( std::string::size_type i = 0; i < s.size(); ++i ) {
	 escape_impl( nrv, s[ i ] );
      }
      return nrv;
   }

   template< int ... Cs > struct escaper;

   template<>
   struct escaper<>
   {
      static std::string result()
      {
	 return std::string();
      }
   };

   template< int C, int ... Cs >
   struct escaper< C, Cs ... >
   {
      static std::string result()
      {
	 return escape( C ) + escaper< Cs ... >::result();
      }
   };

   namespace nocopy_impl
   {
      template< typename >
      class nocopy
      {
      protected:
	 nocopy() {}
	 ~nocopy() {}

      private:
	 nocopy( const nocopy & );
	 void operator= ( const nocopy & );
      };

   } // nocopy_impl

   using namespace nocopy_impl;

   template< typename T > class freer : private nocopy< freer< T > >
   {
   public:
      freer( T * const p )
	    : m_p( p )
      { }

      ~freer()
      {
	 ::free( const_cast< void * >( reinterpret_cast< const volatile void * >( m_p ) ) );
      }

      T * get() const
      {
	 return m_p;
      }

   private:
      T * const m_p;
   };

   inline std::string demangle_impl( const char * const mangled )
   {
      if ( ! mangled ) {
	 return "(null)";
      }
      else if ( ! mangled[ 0 ] ) {
	 return "(empty)";
      }
      int status = -1;
      const freer< const char > demangled( abi::__cxa_demangle( mangled, 0, 0, & status ) );

      if ( ! demangled.get() ) {
	 return mangled;
      }
      else if ( 0 == status ) {
	 return std::string( demangled.get() );
      }
      return mangled;
   }

   template< typename T >
   std::string demangle()
   {
      return demangle_impl( typeid( T ).name() );
   }

   template< typename T >
   std::string demangle( const char *, const char *, const char * )
   {
      return demangle< T >();
   }

   template< typename T1, typename T2, typename ... Ts >
   std::string demangle( const char * a, const char * b, const char * c )
   {
      return a + demangle< T1 >() + b + demangle< T2, Ts ... >( "", b, "" ) + c;
   }

   class file_reader : private nocopy< file_reader >
   {
   public:
      explicit
      file_reader( const std::string & filename )
	    : m_fd( ::open( filename.c_str(), O_RDONLY ) ),
	      m_fn( filename )
      {
	 if ( m_fd < 0 ) {
	    PEGTL_THROW( "unable to open() file " << m_fn << " for reading errno " << errno );
	 }
      }

      ~file_reader()
      {
	 ::close( m_fd );
      }

      size_t size() const
      {
	 struct stat st;

	 errno = 0;
	 if ( ::fstat( m_fd, & st ) < 0 ) {
	    PEGTL_THROW( "unable to fstat() file " << m_fn << " descriptor " << m_fd << " errno " << errno );
	 }
	 return size_t( st.st_size );
      }

      template< typename Container >
      std::string read()
      {
	 Container nrv;
	 nrv.resize( size() );

	 errno = 0;
	 if ( nrv.size() && ( ::read( m_fd, & nrv[ 0 ], nrv.size() ) != int( nrv.size() ) ) ) {
	    PEGTL_THROW( "unable to read() file " << m_fn << " descriptor " << m_fd << " errno " << errno );
	 }
	 return nrv;
      }

   public:
      int internal_fd() const
      {
	 return m_fd;
      }

   private:
      const int m_fd;
      const std::string m_fn;
   };

   class file_mapper : private nocopy< file_mapper >
   {
   public:
      explicit
      file_mapper( const std::string & filename )
	    : m_data( 0 ),
	      m_size( 0 )
      {
	 const file_reader tmp( filename );

	 if ( ! ( m_size = tmp.size() ) ) {
	    return;
	 }
	 errno = 0;

	 if ( intptr_t( m_data = static_cast< const char * >( ::mmap( 0, m_size, PROT_READ, MAP_FILE | MAP_PRIVATE, tmp.internal_fd(), 0 ) ) ) == -1 ) {
	    PEGTL_THROW( "unable to mmap() file " << filename << " errno " << errno );
	 }
      }

      ~file_mapper()
      {
	 ::munmap( const_cast< char * >( m_data ), m_size );
      }

      size_t size() const
      {
	 return m_size;
      }

      const char * data() const
      {
	 return m_data;
      }

      typedef const char * iterator;

      iterator begin() const
      {
	 return m_data;
      }

      iterator end() const
      {
	 return m_data + m_size;
      }

   private:
      const char * m_data;
      size_t m_size;
   };

   inline std::string read_string( const std::string & filename )
   {
      return file_reader( filename ).read< std::string >();
   }

   template< typename Type >
   inline std::vector< Type > read_vector( const std::string & filename )
   {
      return file_reader( filename ).read< std::vector< Type > >();
   }

   template< typename T >
   std::string to_string( const T & t )
   {
      std::ostringstream oss;
      oss << t;
      return oss.str();
   }

   template< class S > S string_to_signed( const std::string & value )
   {
      S nrv = 0;
      S tmp;

      if ( value.empty() ) {
	 PEGTL_THROW( "string-to-integer conversion failed for empty string" );
      }
      unsigned run = 0;

      bool neg = value[ run ] == '-';

      if ( value[ run ] == '+' || value[ run ] == '-' ) {
	 if ( value.size() < 2 ) {
	    PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- sign without number" );
	 }
	 else {
	    ++run;
	 }
      }
      if ( ! isdigit( value[ run ] ) ) {
	    PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- no digits" );
      }
      for ( ; ( run < value.size() ) && ::isdigit( value[ run ] ); ++run )
      {
	 if ( ( tmp = nrv * 10 ) / 10 != nrv ) {
	    PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- integer overflow" );
	 }
	 if ( neg ) {
	    if ( ( nrv = tmp - ( value[ run ] - '0' ) ) + ( value[ run ] - '0' ) != tmp ) {
	       PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- integer overflow" );
	    }
	 }
	 if ( ! neg ) {
	    if ( ( nrv = tmp + ( value[ run ] - '0' ) ) - ( value[ run ] - '0' ) != tmp ) {
	       PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- integer overflow" );
	    }
	 }
      }
      for ( ; run < value.size(); ++run ) {
	 if ( ! ::isspace( value[ run ] ) ) {
	    PEGTL_THROW( "string-to-integer conversion failed for string \"" << value << "\" -- trailing garbage" );
	 }
      }
      return nrv;
   }

} // pegtl

#endif
