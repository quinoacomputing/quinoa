// Copyright (c) 2008 Dr. Colin Hirsch
// Please see license.txt for license.

#include <pegtl.hh>

namespace calculator
{
   // The first program used during development and debugging
   // of the library that uses actions. It evaluates each command
   // line argument as arithmetic expression consisting of
   // - integers with optional sign,
   // - the four basic arithmetic operations,
   // - grouping brackets.
   // For example input "3 * ( -7 + 9)" yields result 6.

   using namespace pegtl;

   // The state.

   // Canonical use of an evaluation stack, here implemented with a std::vector.

   typedef int value_type;
   typedef std::vector< value_type > stack_type;

   // Helper function that is a "value returning pop() operation" that is not
   // in general exception safe, but fine here since the elements are a POD.

   value_type pull( stack_type & s )
   {
      assert( ! s.empty() );
      value_type nrv( s.back() );
      s.pop_back();
      return nrv;
   }

   // The actions.

   // This action converts the matched sub-string to an integer and pushes it on
   // the stack, which must be its only additional state argument.

   // Deriving from action_base<> is necessary since version 0.26; the base class
   // takes care of the pretty-printing for diagnostic messages; this is necessary
   // for all action classes (that do not derive from a rule class).

   struct push_action
	 : action_base< push_action >
   {
      static void apply( const std::string & m, stack_type & s )
      {
	 s.push_back( string_to_signed< value_type >( m ) );
      }
   };

   // Class op_action performs an operation on the two top-most elements of
   // the evaluation stack. This should always be possible in the sense that
   // the grammar must make sure to only apply this action when sufficiently
   // many operands are on the stack.

   template< typename Operation >
   struct op_action
	 : action_base< op_action< Operation > >
   {
      static void apply( const std::string &, stack_type & s )
      {
	 const value_type rhs = pull( s );
	 const value_type lhs = pull( s );
	 s.push_back( Operation()( lhs, rhs ) );
      }
   };

   // Specialisation for division that checks for division by zero.

   template<>
   struct op_action< std::divides< value_type > >
	 : action_base< op_action< std::divides< value_type > > >
   {
      template< typename State >
      static void apply( const std::string &, State & s )
      {
	 const value_type rhs = pull( s );
	 if ( ! rhs ) {
	    PEGTL_THROW( "pegtl: division by zero" );
	 }
	 const value_type lhs = pull( s );
	 s.push_back( std::divides< value_type >()( lhs, rhs ) );
      }
   };

   // Apart for the actions this grammar is identical to the one in parsetree.cc.
   // There, a syntax or expression tree is generated, here the expressions are
   // evaluated on the fly.

   struct read_number
	 : seq< opt< one< '+', '-' > >, plus< digit > > {};

   // This rule uses the rule read_number to match a number in the input
   // and, on success, applies the push_action to the matched sub-string
   // and the state. Here in calculator.cc the state is an instance of
   // stack_type, the evaluation stack on which the number is pushed.

   struct push_rule
	 : pad< ifapply< read_number, push_action >, space > {};

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

   // The operation rules first read two sub-expressions, and then combine their
   // respective results with an arithmetic operation, replacing the two top-most
   // stack elements with the result -- just as every simple stack machine does.

   template< int P, typename O, typename A >
   struct read_op
	 : ifapply< ifmust< calc_pad< P >, O >, op_action< A > > {};

   struct read_mul
	 : read_op< '*', read_atom, std::multiplies< value_type > > {};

   struct read_div
	 : read_op< '/', read_atom, std::divides< value_type > > {};

   struct read_prod
	 : seq< read_atom, star< sor< read_mul, read_div > > > {};

   struct read_add
	 : read_op< '+', read_prod, std::plus< value_type > > {};

   struct read_sub
	 : read_op< '-', read_prod, std::minus< value_type > > {};

   struct read_expr
	 : seq< read_prod, star< sor< read_add, read_sub > > > {};

   struct read_calc
	 : seq< read_expr, space_until_eof > {};

} // calculator

int main( int argc, char ** argv )
{
   for ( int arg = 1; arg < argc; ++arg ) {
      calculator::stack_type stack;
      pegtl::basic_parse_string< calculator::read_calc >( argv[ arg ], stack );
      assert( stack.size() == 1 );  // This can only trigger if the grammar is incorrect.
      std::cerr << "input " << argv[ arg ] << " result " << stack.front() << "\n";
   }
   return 0;
}
