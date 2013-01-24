// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_RULES_ACTION_HH
#define COHI_PEGTL_RULES_ACTION_HH


namespace pegtl
{
   template< typename ... Funcs > struct apply_helper;

   template<>
   struct apply_helper<>
   {
      template< typename ... States >
      static void apply( const std::string &, States && ... )
      { }
   };

   template< typename Func, typename ... Funcs >
   struct apply_helper< Func, Funcs ... >
   {
      template< typename ... States >
      static void apply( const std::string & s, States && ... st )
      {
	 Func::apply( s, std::forward< States >( st ) ... );
	 apply_helper< Funcs ... >::apply( s, std::forward< States >( st ) ... );
      }
   };

   template< typename Func, typename ... Funcs >
   struct apply
   {
      typedef apply key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< apply, Func, Funcs ... >( st, "$", "", " $", "", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input &, Debug &, States && ... st )
      {
	 apply_helper< Func, Funcs ... >::apply( "", std::forward< States >( st ) ... );
	 return true;
      }
   };

   template< typename Rule, typename Func, typename ... Funcs >
   struct ifapply
   {
      typedef ifapply key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare2< ifapply, Rule, Func, Funcs ... >( st, "( ", " | $", " $", " )" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< false > p( in );

	 if ( Rule::template match< Must >( in, de, std::forward< States >( st ) ... ) ) {
	    apply_helper< Func, Funcs ... >::apply( std::string( p.here(), in.here() ), std::forward< States >( st ) ... );
	    return p( true );
	 }
	 return p( false );
      }
   };

   template< typename Func >
   struct action_base
   {
      typedef Func key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< Func >( demangle< Func >(), true );
      }
   };

   struct nop
   {
      typedef nop key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< nop >( "nop", true );
      }

      template< typename ... States >
      static void apply( const std::string &, States && ... )
      { }
   };

   template< typename Rule >
   struct ifapply< Rule, nop >
	 : Rule {};

   template< unsigned N, typename ... Funcs >
   struct nth
   {
      typedef nth key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = to_string( N );
	 prepare1< nth, Funcs ... >( st, ( n + ":" ).c_str(), "", ( " $" + n + ":" ).c_str(), "", "" );
      }

      template< typename State, typename ... States >
      static void apply( const std::string & s, State &&, States && ... st )
      {
	 nth< N - 1, Funcs ... >::template apply< States ... >( s, std::forward< States >( st ) ... );
      }
   };

   template< typename ... Funcs >
   struct nth< 0, Funcs ... >
   {
      typedef nth key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< nth, Funcs ... >( st, "0:", "", " $0:", "", "" );
      }

      template< typename State, typename ... States >
      static void apply( const std::string & s, State && st, States && ... )
      {
	 apply_helper< Funcs ... >::apply( s, std::forward< State >( st ) );
      }
   };

   struct insert
   {
      typedef insert key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< nop >( "insert", true );
      }

      template< typename Container >
      static void apply( const std::string & s, Container & c )
      {
	 c.insert( c.end(), s );
      }
   };

   struct push_back
   {
      typedef push_back key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< nop >( "push_back", true );
      }

      template< typename Container >
      static void apply( const std::string & s, Container & c )
      {
	 c.push_back( s );
      }
   };

   // Class capture is both an action and a rule.

   // As action, it stores the matched string as value in the
   // 'capture_map &' state, using the template argument as key.

   // As rule, it retrieves a string from the 'capture_map' state,
   // again using the template argument as key, and attempts to
   // match the input against that string.

   template< unsigned Key >
   struct capture
   {
      typedef capture key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< capture >( '\\' + to_string( Key ), true );
      }

      template< typename Map >
      static void apply( const std::string & s, Map & m )
      {
	 m[ Key ] = s;
      }

      template< bool Must, typename Input, typename Debug, typename Map >
      static bool match( Input & in, Debug &, const Map & m )
      {
	 typename Input::template marker< Must > p( in );

	 const typename Map::const_iterator i = m.find( Key );

	 if ( i == m.end() ) {
	    return p( false );  // Or assume empty string? No, too dangerous wrt. infinite recursions.
	 }
	 const std::string & s = i->second;

	 for ( size_t j = 0; j < s.size(); ++j ) {
	    character< Input > c( in );
	    if ( ! c( c == s[ j ] ) ) {
	       return p( false );
	    }
	 }
	 return p( true );
      }
   };

   // Typedef for a map that can be used together with class capture.

   typedef std::map< unsigned, std::string > capture_map;

} // pegtl

#endif
