// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace example
{
   // This file shows how to create a syntax tree for simple
   // arithmetic expressions similar to what calculator.cc
   // handles. The basic approach is very stack-machine like,
   // again following the lead from calculator.cc, but this
   // time to create a data structure, rather than immediately
   // calculate a result.

   using namespace pegtl;

   // The node class for the syntax tree. It is used both as
   // inner node, in which case the sons data member is used,
   // and as leaf node, in which case sons is empty and the
   // matched data member contains the matched sub-string.

   struct tree_node
   {
      const unsigned id;
      const std::string matched;
      const std::vector< std::shared_ptr< tree_node > > sons;

      tree_node( const int id, const std::string & matched )
	    : id( id ),
	      matched( matched )
      { }

      template< typename Iterator >
      tree_node( const int id, const Iterator & begin, const Iterator & end )
	    : id( id ),
	      sons( begin, end )
      { }
   };

   // This stack is used as state in the parser; during the parse,
   // the actions use the stack for intermediate results on which
   // they operate in a stack-machine style. At the end of a successful
   // parse the stack contains the resulting syntax tree.

   typedef std::vector< std::shared_ptr< tree_node > > tree_stack;

   // Action class for leaf nodes; it simply pushes a new with the
   // matched string onto the stack.

   template< unsigned Id >
   struct push_node
	 : action_base< push_node< Id > >
   {
      static void apply( const std::string & s, tree_stack & t )
      {
	 t.push_back( std::make_shared< tree_node >( Id, s ) );
      }
   };

   // Action class for inner nodes; it assumes that there are Size
   // elements present on the stack, which are then used as sons
   // of a newly created inner node; the newly created node replaces
   // the nodes that have become its sons on the stack.

   template< unsigned Id, unsigned Size >
   struct make_node
	 : action_base< make_node< Id, Size > >
   {
      static void apply( const std::string &, tree_stack & t )
      {
	 assert( t.size() >= Size );
	 const std::shared_ptr< tree_node > x = std::make_shared< tree_node >( Id, t.end() - Size, t.end() );
	 t.erase( t.end() - Size, t.end() );
	 t.push_back( x );
      }
   };

   // Output operator for syntax trees, for debugging purposes...

   std::ostream & operator<< ( std::ostream & o, const std::shared_ptr< tree_node > & t )
   {
      if ( ! t ) {
	 o << "( )";
	 return o;
      }
      o << "id=" << t->id
	<< " matched='" << t->matched
	<< "' sons=" << t->sons.size()
	<< " (";
      for ( std::vector< std::shared_ptr< tree_node > >::const_iterator i = t->sons.begin(); i != t->sons.end(); ++i ) {
	 o << " " << ( * i );
      }
      o << " )";
      return o;
   }

   // Apart for the actions this grammar is identical to the one in calculator.cc.
   // There, the expressions are evaluated on the fly, here a syntax or expression
   // tree is generated.

   struct read_number
         : seq< opt< one< '+', '-' > >, plus< digit > > {};

   // This rule uses the rule read_number to match a number in the input
   // and, on success, applies the push_node action to the matched sub-
   // string which, like in calculator.cc, pushes a new element onto
   // the stack; here as leaf node, there as integer.

   struct push_rule
         : pad< ifapply< read_number, push_node< 0 > >, space > {};

   template< int C >
   struct calc_pad
         : pad_one< C, space > {};

   struct read_open
         : calc_pad< '(' > {};

   struct read_close
         : calc_pad< ')' > {};

   struct read_expr;

   struct read_atom
         : sor< push_rule, seq< read_open, read_expr, read_close > > {};

   // The operation rules operate on the stack; here the two elements
   // popped off the stack are replaced by a new inner node that points
   // to the two popped elements, in calculator.cc the operation was
   // actually performed and the two popped elements replaced by the
   // result of the operation.

   template< unsigned P, typename O >
   struct read_op
	 : ifmust< calc_pad< P >, O, apply< make_node< P, 2 > > > {};

   struct read_mul
	 : read_op< '*', read_atom > {};

   struct read_div
	 : read_op< '/', read_atom > {};

   struct read_prod
         : seq< read_atom, star< sor< read_mul, read_div > > > {};

   struct read_add
	 : read_op< '+', read_prod > {};

   struct read_sub
	 : read_op< '-', read_prod > {};

   struct read_expr
         : seq< read_prod, star< sor< read_add, read_sub > > > {};

   struct grammar
         : until< space_until_eof, read_expr > {};

} // example

int main( int argc, char ** argv )
{
   for ( int i = 1; i < argc; ++i ) {

      // To use the grammar, an empty stack must be passed as state.

      example::tree_stack s;
      pegtl::basic_parse_string< example::grammar >( argv[ i ], s );
      std::cout << "parsed " << s.size() << " top-level expressions\n";

      // The for-loop is necessary because the grammar accepts arbitrary
      // many full expression before the eof, each of which will leave one
      // pointer to the root of the corresponding syntax tree in the stack.

      for ( size_t j = 0; j < s.size(); ++j ) {
	 std::cout << s[ j ] << "\n";
      }
   }
   return 0;
}
