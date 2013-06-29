// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_RULES_SPECIAL_HH
#define COHI_PEGTL_RULES_SPECIAL_HH

// #define PEGTL_IMPURE_OPTIMISATIONS

namespace pegtl
{
   template<>
   struct at< success >
	 : success {};

   template<>
   struct at< failure >
	 : failure {};

   template<>
   struct not_at< success >
	 : failure {};

   template<>
   struct not_at< failure >
	 : success {};

   template<>
   struct opt< success >
	 : success {};

   template<>
   struct opt< failure >
	 : success {};

   template<>
   struct plus< success >
	 : success {};

   template<>
   struct plus< failure >
	 : failure {};

   template<>
   struct star< success >
	 : success {};

   template<>
   struct star< failure >
	 : success {};

   template<>
   struct at< eof >
	 : eof {};

   template< typename Rule >
   struct at< at< Rule > >
	 : at< Rule > {};

   template< typename Rule >
   struct not_at< not_at< Rule > >
	 : at< Rule > {};

   template< typename Rule >
   struct not_at< at< Rule > >
	 : not_at< Rule > {};

   template< typename Rule >
   struct at< not_at< Rule > >
	 : not_at< Rule > {};

#if defined( PEGTL_IMPURE_OPTIMISATIONS )

   template< typename Rule >
   struct opt< at< Rule > >
	 : success {};

   template< typename Rule >
   struct opt< not_at< Rule > >
	 : success {};

   template< typename Rule >
   struct at< opt< Rule > >
	 : success {};

   template< typename Rule >
   struct not_at< opt< Rule > >
	 : success {};

   template< typename Rule >
   struct star< at< Rule > >
	 : success {};

   template< typename Rule >
   struct star< not_at< Rule > >
	 : success {};

   template< typename Rule >
   struct at< star< Rule > >
	 : success {};

   template< typename Rule >
   struct not_at< star< Rule > >
	 : success {};

   template< typename Rule >
   struct plus< at< Rule > >
	 : at< Rule > {};

   template< typename Rule >
   struct plus< not_at< Rule > >
	 : not_at< Rule > {};

   template< typename Rule >
   struct at< plus< Rule > >
	 : at< Rule > {};

   template< typename Rule >
   struct not_at< plus< Rule > >
	 : not_at< Rule > {};

   template< typename Rule >
   struct opt< opt< Rule > >
	 : opt< Rule > {};

   template< typename Rule >
   struct plus< plus< Rule > >
	 : plus< Rule > {};

   template< typename Rule >
   struct opt< plus< Rule > >
	 : star< Rule > {};

   template< typename Rule >
   struct opt< star< Rule > >
	 : star< Rule > {};

   template< typename Rule >
   struct star< plus< Rule > >
	 : star< Rule > {};

#endif // PEGTL_IMPURE_OPTIMISATIONS

   template< typename Rule >
   struct star< opt< Rule > >
   {
      static_assert( !sizeof( Rule ), "pegtl: illegal expression Rule?* (allows iteration without progress = infinite loop)" );
   };

   template< typename Rule >
   struct plus< opt< Rule > >
   {
      static_assert( !sizeof( Rule ), "pegtl: illegal expression Rule?+ (allows iteration without progress = infinite loop)" );
   };

   template< typename Rule >
   struct star< star< Rule > >
   {
      static_assert( !sizeof( Rule ), "pegtl: illegal expression Rule** (allows iteration without progress = infinite loop)" ); 
   };

   template< typename Rule >
   struct plus< star< Rule > >
   {
      static_assert( !sizeof( Rule ), "pegtl: illegal expression Rule*+ (allows iteration without progress = infinite loop)" );
   };

   template< typename Rule, typename PadL, typename PadR = PadL >
   struct pad
         : rule_base< pad< Rule, PadL, PadR >, seq< star< PadL >, Rule, star< PadR > > > {};

   template< typename Rule, typename PadL >
   struct padl
         : rule_base< padl< Rule, PadL >, seq< star< PadL >, Rule > > {};

   template< typename Rule, typename PadR >
   struct padr
         : rule_base< padr< Rule, PadR >, seq< Rule, star< PadR > > > {};

   template< typename Rule, typename Glue, typename Action = nop >
   struct list
         : rule_base< list< Rule, Glue, Action >, seq< Rule, star< ifmust< Glue, ifapply< Rule, Action > > > > > {};

   template< typename Begin, typename Body, typename End = Begin, typename Action = nop >
   struct enclose
	 : rule_base< enclose< Begin, Body, End, Action >, ifmust< Begin, ifapply< until< at< End >, Body >, Action >, End > > {};

} // pegtl

#endif
