// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_RULES_BASIC_HH
#define COHI_PEGTL_RULES_BASIC_HH


namespace pegtl
{
   template< bool B >
   struct bool_rule
   {
      typedef bool_rule key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 const std::string n = B ? "success" : "failure";
	 st.template update< bool_rule >( n, true );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input &, Debug &, States && ... )
      {
	 return B;
      }
   };

   struct success
	 : bool_rule< true > {};

   struct failure
	 : bool_rule< false > {};

   template< typename ... > struct seq;

   template<>
   struct seq<>
	 : success
   {
      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool i_match( Input &, Debug &, States && ... )
      {
	 return true;
      }
   };

   template< typename Rule, typename ... Rules >
   struct seq< Rule, Rules ... >
   {
      typedef seq key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< seq, Rule, Rules ... >( st, "", "( ", " ", " )", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< Must || ! sizeof ... ( Rules ) > h( in );
	 return h( i_match< Must >( in, de, std::forward< States >( st ) ... ) );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool i_match( Input & in, Debug & de, States && ... st )
      {
	 return de.template match< Must, Rule >( in, std::forward< States >( st ) ... ) && pegtl::seq< Rules ... >::template i_match< Must >( in, de, std::forward< States >( st ) ... );
      }
   };

   template< typename ... Rules >
   struct opt
   {
      typedef opt key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< opt, Rules ... >( st, "", "( ", " ", " )", "?" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 return in.eof() || seq< Rules ... >::template match< false >( in, de, std::forward< States >( st ) ... ) || true;
      }
   };

   template< typename ... Rules >
   struct star
   {
      typedef star key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< star, Rules ... >( st, "", "( ", " ", " )", "*" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 while ( ( ! in.eof() ) && seq< Rules ... >::template match< false >( in, de, std::forward< States >( st ) ... ) ) {}
	 return true;
      }
   };

   template< typename ... Rules >
   struct plus
   {
      typedef plus key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< plus, Rules ... >( st, "", "( ", " ", " )", "+" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 return seq< Rules ... >::template match< Must >( in, de, std::forward< States >( st ) ... ) && star< Rules ... >::template match< Must >( in, de, std::forward< States >( st ) ... );
      }
   };

   template< typename... > struct sor;

   template<>
   struct sor<>
	 : failure {};

   template< typename Rule, typename... Rules >
   struct sor< Rule, Rules ... >
   {
      typedef sor key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< sor, Rule, Rules ... >( st, "", "( ", " / ", " )", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 return de.template match< false, Rule >( in, std::forward< States >( st ) ... ) || pegtl::sor< Rules ... >::template match< false >( in, de, std::forward< States >( st ) ... );
	 // return de.template match< Must && ! sizeof ... ( Rules ), Rule >( in, std::forward< States >( st ) ... ) || pegtl::sor< Rules ... >::template match< Must >( in, de, std::forward< States >( st ) ... );
      }
   };

   template< typename Rule >
   struct at
   {
      typedef at key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< at, Rule >( st, "&", "( ", " ", " )", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< false > p( in );
	 return de.template match< Must, Rule >( in, std::forward< States >( st ) ... );
      }
   };

   template< typename Rule >
   struct not_at
   {
      typedef not_at key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 prepare1< not_at, Rule >( st, "!", "( ", " ", " )", "" );
      }

      template< bool Must, typename Input, typename Debug, typename ... States >
      static bool match( Input & in, Debug & de, States && ... st )
      {
	 typename Input::template marker< false > p( in );
	 return ! de.template match< false, Rule >( in, std::forward< States >( st ) ... );
      }
   };

} // pegtl

#endif
