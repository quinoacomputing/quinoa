// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_RULES_EXTENDED_HH
#define COHI_PEGTL_RULES_EXTENDED_HH


namespace pegtl
{
   template< unsigned N, typename ... Rules >
   struct rep
   {
      typedef rep key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< rep, Rules ... >( st, "", "( ", " ", " )", ( "{" + to_string( N ) + "}" ).c_str() );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must > p( in );

	 for ( unsigned i = 0; i < N; ++i ) {
	    if ( ! seq< Rules ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... ) ) {
	       return p( false );
	    }
	 }
	 return p( true );
      }
   };

   template< typename ... Rules >
   struct must
   {
      typedef must key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< must, Rules ... >( st, "!", "( ", " ", " )", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 return seq< Rules ... >::template match< true >( in, de, std::forward< States >( st ) ... );
      }
   };

   struct eof
   {
      typedef eof key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template update< eof >( "&eof", true );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug &, States && ... )
      {
	 return in.eof();
      }
   };

   template< typename Cond, typename ... Rules > struct until;

   template< typename Cond, typename Rule, typename ... Rules >
   struct until< Cond, Rule, Rules ... >
   {
      typedef until key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template insert< Cond >();
	 const std::string n = st.template name< Cond >();
	 prepare1< until, Rule, Rules ... >( st, "", "( ", " ", " )", ( "{ " + n + " }" ).c_str());
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must > p( in );

	 while ( ! de.template match< false, Cond >( in, std::forward< States >( st ) ... ) ) {
	    if ( in.eof() || ( ! seq< Rule, Rules ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... ) ) ) {
	       return p( false );
	    }
	 }
	 return p( true );
      }
   };

   template< typename Cond >
   struct until< Cond >
   {
      typedef until key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< until, Cond >( st, ".{ ", "", "", "", " }" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
         typename Input::template marker< Must > p( in );

         while ( ! de.template match< false, Cond >( in, std::forward< States >( st ) ... ) ) {
            if ( in.eof() ) {
               return p( false );
            }
            in.bump();
         }
         return p( true );
      }
   };

   template< typename Rule, typename Impl >
   struct rule_base
	 : Impl
   {
      template< typename Print >
      static void prepare( Print & st )
      {
	 Impl::prepare( st );
	 const std::string y = demangle< Rule >();
	 const std::string m = st.template name< Rule >();
	 const std::string e = st.template expr< Impl >();
	 st.template update< Rule >( e, m == y );
      }
   };

   template< bool Must, typename Cond, typename ... Thens >
   struct cond2impl
   {
      typedef cond2impl key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare2< cond2impl, Cond, Thens ... >( st, "( ", Must ? " => " : " -> ", " ", " )" );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must > p( in );

	 if ( de.template match< false, Cond >( in, std::forward< States >( st ) ... ) ) {
	    return p( seq< Thens ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... ) );
	 } else {
	    return p( ! Must );
	 }
      }
   };

   template< typename Cond, typename ... Thens >
   struct ifmust
	 : rule_base< ifmust< Cond, Thens ... >, cond2impl< true, Cond, Thens ... > > {};

   template< typename Cond, typename ... Thens >
   struct ifthen
	 : rule_base< ifthen< Cond, Thens... >, cond2impl< false, Cond, Thens ... > > {};

   template< bool Must, typename Cond, typename Then, typename Else >
   struct cond3impl
   {
      typedef cond3impl key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
         prepare2< cond3impl, Cond, Then, Else >( st, "( ", Must ? " => " : " -> ", " / ", " )" );
      }

      template< bool, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must > p( in );

	 if ( de.template match< false, Cond >( in, std::forward< States >( st ) ... ) ) {
	    return p( de.template match< Must, Then >( in, std::forward< States >( st ) ... ) );
	 }
	 else {
	    return p( de.template match< Must, Else >( in, std::forward< States >( st ) ... ) );
	 }
      }
   };

   template< typename Cond, typename Then, typename Else >
   struct ifmustelse
	 : rule_base< ifmustelse< Cond, Then, Else >, cond3impl< true, Cond, Then, Else > > {};

   template< typename Cond, typename Then, typename Else >
   struct ifthenelse
	 : rule_base< ifthenelse< Cond, Then, Else >, cond3impl< false, Cond, Then, Else > > {};

   template< unsigned N, unsigned M, typename ... Rules >
   struct rep2
   {
      typedef rep2 key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< rep2, Rules ... >( st, "", "( ", " ", " )", ( "{" + to_string( N ) + ',' + to_string( M ) + "}" ).c_str() );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 static_assert( N <= M, "pegtl: illegal expression rep< R, N, M > where N is greater than M" );

	 typename Input::template marker< Must > p( in );

	 for ( unsigned i = 0; i < N; ++i ) {
	    if ( ! seq< Rules ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... ) ) {
	       return p( false );
	    }
	 }
	 for ( unsigned i = N; i < M; ++i ) {
	    if ( ! seq< Rules ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... ) ) {
	       return p( true );
	    }
	 }
	 return p( not_at< seq< Rules ... > >::template match< Must >( in, de, std::forward< States >( st ) ... ) );
      }
   };

} // pegtl

#endif
