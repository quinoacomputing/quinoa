// Copyright (c) 2008 by Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace example
{
   using namespace pegtl;

   // First example on how to parse a cons-style list. It uses an
   // intrusive style where the main action for building the lists
   // is not an action in the usual pegtl terminology sense. Rather
   // it is a (re-)implementation of a parsing rule that contains
   // the cons-list building code embedded at the appropriate places.

   // + All operations O(1)
   // - C++ stack-usage proportional to the nesting depth of the cons-cells
   // + Only a single action/rule to handle the list
   // - Intrusive with duplicated parsing rules with embedded tree-building
   // + Very simple state consists of a single node reference
   // + Cons-cells can have car and cdr declared const

   // Simple classes for the representation of cons-style
   // lists just enough to show the principles...

   // Base class for all other classes.

   struct node_base
	 : private nocopy< node_base >
   {
      virtual ~node_base() {}

      virtual void v_print( std::ostream & ) const = 0;
   };

   struct cons_node;

   std::ostream & operator<< ( std::ostream & o, const std::shared_ptr< cons_node > & p );
   std::ostream & operator<< ( std::ostream & o, const std::shared_ptr< node_base > & p )
   {
      if ( p ) {
	 p->v_print( o );
      }
      else {
	 o << "nil";
      }
      return o;
   }

   // Class for cons-style list; the cdr is a pointer to cons_node rather than node_base
   // because it simplifies the examples -- that only generate proper lists, i.e. where
   // the cdr is either another cons-cell or nil/null.

   struct cons_node
	 : public node_base
   {
      cons_node( const std::shared_ptr< node_base > & a,
		 const std::shared_ptr< node_base > & d )
	    : m_car( a ),
	      m_cdr( d )
      { }

      const std::shared_ptr< node_base > m_car;
      const std::shared_ptr< node_base > m_cdr;

      void v_print( std::ostream & o ) const
      {
	 o << "( " << m_car;

	 // Hackish, the grammar implicitly ensures that the cast always succeeds wherever necassary...

	 for ( std::shared_ptr< cons_node > t = std::dynamic_pointer_cast< cons_node >( m_cdr ); t; t = std::dynamic_pointer_cast< cons_node >( t->m_cdr ) ) {
	    o << " " << t->m_car;
	 }
	 o << " )";
      }
   };

   std::ostream & operator<< ( std::ostream & o, const std::shared_ptr< cons_node > & p )
   {
      if ( p ) {
	 p->v_print( o );
      }
      else {
	 o << "nil";
      }
      return o;
   }

   // Class for an atomic object; atomic in the sense that it doesn't contain
   // references to other objects; here simply a string.

   struct token_node
	 : public node_base
   {
      explicit
      token_node( const std::string & value )
	    : m_value( value )
      { }

      const std::string m_value;

      void v_print( std::ostream & o ) const
      {
	 o << m_value;
      }
   };

   // The action classes.

   struct make_null
	 : action_base< make_null >
   {
      static void apply( const std::string &, std::shared_ptr< node_base > & result )
      {
	 result.reset();
      }
   };

   // This class is an action that, in the grammar below, is associated with the rule for
   // a token; it creates a new token and "returns" it in the passed-in "result"; this is
   // in line with the whole architecture: The caller of a rule must provide the return
   // value, and pass it as argument. It ignores that stack that is only used for lists.

   struct make_token
	 : action_base< make_token >
   {
      static void apply( const std::string & matched, std::shared_ptr< node_base > & result )
      {
	 result = std::make_shared< token_node >( matched );
      }
   };

   // This class is a rule that contains the necessary operations to create a cons-cell.
   // It is quite straightforward, first it parses the head, then the tail, and then it
   // creates a cons-cell that collects the two. Disregarding the cons-cell-building, it
   // is of course equivalent to seq< Head, Tail >, and it is this redundancy that makes
   // the approach ugly.

   template< typename Head, typename Tail >
   struct make_cons
   {
      typedef seq< Head, Tail > key_type;

      template< typename Print >
      static void prepare( Print & st )
      {
	 st.template insert< Head, Tail >();
      }

      template< bool Must, typename Input, typename Debug >
      static bool match( Input & in, Debug & de, std::shared_ptr< node_base > & result )
      {
	 typename Input::template marker< Must > p( in );

	 std::shared_ptr< node_base > car;
	 std::shared_ptr< node_base > cdr;

	 if ( ! de.template match< Must, Head >( in, car ) ) {
	    return p( false );
	 }
	 if ( ! de.template match< Must, Tail >( in, cdr ) ) {
	    return p( false );
	 }
	 result = std::make_shared< cons_node >( car, cdr );
	 return p( true );
      }
   };

   struct comment
	 : ifmust< one< ';' >, until< eol > > {};

   struct separator
	 : sor< comment, space > {};

   struct read_atom
	 : pad< ifapply< plus< digit >, make_token >, separator > {};

   struct read_expr;
   struct read_tail;

   struct read_cons
	 : make_cons< read_expr, read_tail > {};

   struct read_tail
	 : sor< pad_one< ')', separator >, read_cons > {};

   struct read_list
	 : ifmust< pad_one< '(', separator >, sor< ifapply< one< ')' >, make_null >, read_cons > > {};

   struct read_expr
	 : sor< read_list, read_atom > {};

   struct read_file
	 : seq< read_expr, until< eof, separator > > {};

} // example

int main( int argc, char ** argv )
{
   for ( int arg = 1; arg < argc; ++arg ) {
      std::shared_ptr< example::node_base > result;
      pegtl::basic_parse_string< example::read_file >( argv[ arg ], result );
      PEGTL_PRINT( "input: " << argv[ arg ] << "\nresult: " << result );
   }
   return 0;
}
