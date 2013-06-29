// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace example
{
   using namespace pegtl;

   // Second example on how to parse a cons-style list. It uses the
   // stack-machine style that turning out to be a rather useful way
   // of approaching the "how to attach actions to a grammar" question.

   // + All operations O(1)
   // + C++ stack-usage proportional to the nesting depth of the lists
   // - Requires three distinct actions to handle the list
   // + Non-intrusive with actions attached via apply and ifapply
   // - More complex state consists of a node reference and a stack
   // - Cons-cells must be mutable (to prevent a const_cast)

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
      explicit
      cons_node( const std::shared_ptr< node_base > & a = std::shared_ptr< node_base >(),
		 const std::shared_ptr< cons_node > & d = std::shared_ptr< cons_node >() )
	    : m_car( a ),
	      m_cdr( d )
      { }

      std::shared_ptr< node_base > m_car;
      std::shared_ptr< cons_node > m_cdr;

      void v_print( std::ostream & o ) const
      {
	 o << "( " << m_car;

	 for ( std::shared_ptr< cons_node > t = m_cdr; t; t = t->m_cdr ) {
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

   // Stack type that is used to keep track of the nesting of lists.
   // The .first entry is the head of the list, i.e. the result.
   // The .second entry is the tail of the list and is updated while the list grows.
   // It exists only to make appending a new element to the list O(1), rather than
   // having to walk the list from beginning to end every time...

   typedef std::vector< std::pair< std::shared_ptr< cons_node >, std::shared_ptr< cons_node > > > t_stack;

   // This class is an action that, in the grammar below, is associated with the rule for
   // a token; it creates a new token and "returns" it in the passed-in "result"; this is
   // in line with the whole architecture: The caller of a rule must provide the return
   // value, and pass it as argument. It ignores that stack that is only used for lists.

   struct make_token
	 : action_base< make_token >
   {
      static void apply( const std::string & matched, std::shared_ptr< node_base > & result, t_stack & )
      {
	 std::shared_ptr< node_base >( new token_node( matched ) ).swap( result );
      }
   };

   // The first element of a list pushes a new element onto the stack.

   struct list_first
	 : action_base< list_first >
   {
      static void apply( const std::string &, std::shared_ptr< node_base > & result, t_stack & stack )
      {
	 const std::shared_ptr< cons_node > tmp = std::make_shared< cons_node >( result );
	 stack.push_back( std::make_pair( tmp, tmp ) );
	 result.reset();
      }
   };

   // Consecutive elements append to the list.

   struct list_next
	 : action_base< list_next >
   {
      static void apply( const std::string &, std::shared_ptr< node_base > & result, t_stack & stack )
      {
	 stack.back().second->m_cdr = std::make_shared< cons_node >( result );
	 stack.back().second = stack.back().second->m_cdr;
	 result.reset();
      }
   };

   // Finally the stack entry for the list is removed, and the list head placed as result.

   struct list_end
	 : action_base< list_end >
   {
      static void apply( const std::string &, std::shared_ptr< node_base > & result, t_stack & stack )
      {
	 result = stack.back().first;
	 stack.pop_back();
      }
   };

   struct comment
	 : ifmust< one< ';' >, until< eol > > {};

   struct separator
	 : sor< comment, space > {};

   struct read_atom
	 : pad< ifapply< plus< digit >, make_token >, separator > {};

   struct read_expr;

   struct read_tail
	 : until< pad_one< ')', separator >, read_expr, apply< list_next > > {};

   struct read_head
	 : seq< read_expr, apply< list_first >, read_tail, apply< list_end > > {};

   struct read_list
	 : ifmust< pad_one< '(', separator >, sor< one< ')' >, read_head > > {};

   struct read_expr
	 : sor< read_list, read_atom > {};

   struct read_file
	 : seq< read_expr, until< eof, separator > > {};

} // example

int main( int argc, char ** argv )
{
   for ( int arg = 1; arg < argc; ++arg ) {
      example::t_stack stack;
      std::shared_ptr< example::node_base > result;
      pegtl::basic_parse_string< example::read_file >( argv[ arg ], result, stack );
      PEGTL_PRINT( "input: " << argv[ arg ] << "\nresult: " << result );
   }
   return 0;
}
