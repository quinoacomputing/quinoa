// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#ifndef COHI_PEGTL_HH
#error "Please #include only pegtl.hh (rather than individual pegtl_*.hh files)."
#endif

#ifndef COHI_PEGTL_DEBUG_DUMMY_HH
#define COHI_PEGTL_DEBUG_DUMMY_HH


namespace pegtl
{
   struct debug_base : private nocopy< debug_base >
   {
   public:
      virtual ~debug_base()
      { }

      virtual void log( const std::string & message )
      {
	 std::cerr << "pegtl: " << message << std::endl;
      }

   protected:
      debug_base()
      { }
   };

   struct dummy_debug : public debug_base
   {
      template< bool Must, typename Rule, typename Input, typename ... States >
      bool match( Input & in, States && ... st )
      {
	 if ( Rule::template match< Must >( in, * this, std::forward< States >( st ) ... ) ) {
	    return true;
	 }
	 else if ( Must ) {
	    PEGTL_THROW( "required rule " << demangle< Rule >() << " failed" );
	 }
	 return false;
      }
   };

} // pegtl

#endif
