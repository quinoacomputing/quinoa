// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_DEBUG_PRINT_HH
#define COHI_PEGTL_DEBUG_PRINT_HH

namespace pegtl
{
   class printer
   {
   private:
      typedef void ( printer::*inserter ) ( bool );

      struct value_type
      {
	 value_type()
	 { }

	 value_type( const std::string & name, const std::string & expr )
	       : m_name( name ),
		 m_expr( expr )
	 { }

	 std::string m_name;
	 std::string m_expr;
      };

      typedef std::string key_type;
      typedef std::map< key_type, value_type > map_type;

   public:
      template< typename TopRule >
      printer( const tag< TopRule > & )
               : m_top_inserter( & printer::insert< TopRule > )
      { }

      template< typename Rule >
      static const char* key()
      {
	 return typeid( typename Rule::key_type ).name();
      }

      template< typename Rule >
      static value_type value()
      {
	 const std::string& d = demangle< Rule >();
	 return value_type( d, d );
      }

      template< typename Rule >
      std::string rule()
      {
         return rule( find< Rule >() );
      }

      template< typename Rule >
      const std::string & name()
      {
	 return find< Rule >().m_name;
      }

      template< typename Rule >
      const std::string & expr()
      {
	 return find< Rule >().m_expr;
      }

      template< typename Rule >
      void insert( const bool force = false )
      {
	 if ( m_rules.insert( map_type::value_type( key< Rule >(), value< Rule >() ) ).second || force ) {
	    Rule::prepare( *this );
	 }
      }

      template< typename Rule1, typename Rule2, typename... Rules >
      void insert( const bool force = false )
      {
	 insert< Rule1 >( force );
	 insert< Rule2, Rules... >( force );
      }

      template< typename Rule >
      void update( const std::string & expr, const bool both = false )
      {
         update( key< Rule >(), expr, both );
      }

      void print_rules()
      {
         ensure_insert();

	 for ( map_type::const_iterator i = m_rules.begin(); i != m_rules.end(); ++i ) {
	    if ( i->second.m_name == i->second.m_expr ) {
	       PEGTL_PRINT( "RULE1 " << i->second.m_name );
	    }
	    else {
	       PEGTL_PRINT( "RULE2 " << i->second.m_name << " := " << i->second.m_expr );
	    }
	 }
      }

   private:
      map_type m_rules;

      const value_type m_default_value;

      const inserter m_top_inserter;

      void ensure_insert()
      {
	 if ( m_rules.empty() ) {
	    ( this->*m_top_inserter )( true );
	 }
      }

      template< typename Rule >
      const value_type & find()
      {
         return find( key< Rule >() );
      }

      const value_type & find( const char* key )
      {
         ensure_insert();

         const map_type::const_iterator i = m_rules.find( key );
	 if( i == m_rules.end() ) {
            return m_default_value;
	 }
	 return i->second;
      }

      std::string rule( const value_type& e )
      {
	 if ( e.m_name.empty() || ( e.m_name == e.m_expr ) ) {
	    return e.m_name;
	 }
	 else {
	    return e.m_name + " := " + e.m_expr;
	 }
      }

      void update( const char* key, const std::string & expr, const bool both )
      {
	 if ( both ) {
            m_rules[ key ] = value_type( expr, expr );
	 }
	 else {
	    m_rules[ key ].m_expr = expr;
	 }
      }
   };

   template< typename Rule, typename Print >
   inline const std::string & names( Print & st, const char *, const char *, const char * )
   {
      return st.template name< Rule >();
   }

   template< typename Rule1, typename Rule2, typename ... Rules, typename Print >
   inline std::string names( Print & st, const char * a, const char * b, const char * c )
   {
     std::string nrv = a;
     nrv += names< Rule1 >( st, "", "", "" );
     nrv += b;
     nrv += names< Rule2, Rules ... >( st, "", b, "" );
     nrv += c;
     return nrv;
   }

   template< typename Master, typename Rule, typename ... Rules, typename Print >
   inline void prepare1( Print & st, const char * a, const char * b, const char * c, const char * d, const char * e )
   {
      st.template insert< Rule, Rules ... >();
      std::string s = a;
      s += names< Rule, Rules ... >( st, b, c, d );
      s += e;
      const bool u = demangle< Master >() == st.template name< Master >();
      st.template update< Master >( s , u );
   }

   template< typename Master, typename Head, typename Rule, typename ... Rules, typename Print >
   inline void prepare2( Print & st, const char * a, const char * b, const char * c, const char * d )
   {
      st.template insert< Head, Rule, Rules ... >();
      std::string s = a;
      s += st.template name< Head >();
      s += b;
      s += names< Rule, Rules ... >( st, "", c, "" );
      s += d;
      const bool u = demangle< Master >() == st.template name< Master >();
      st.template update< Master >( s, u );
   }

   template< typename TopRule >
   inline void print_rules()
   {
      printer( tag< TopRule >() ).print_rules();
   }

} // pegtl

#endif
