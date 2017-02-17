/**   ------------------------------------------------------------
 *    Copyright 2003 - 2011 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

/*
  Loosely based on:  (Less and less so each checkin, practically none)
  File: Eval.c
  Auth: Brian Allen Vanderburg II
  Date: Wednesday, April 30, 2003
  Desc: Evaluation routines for the Eval library
*/


#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <cstdlib>
#include <cctype>
#include <cmath>

#include <time.h>

#include <Sacado.hpp>

#include <stk_expreval/ExprEvalFAD.hpp>
#include <stk_expreval/Lexer.hpp>

namespace stk_classic {
namespace expreval {
namespace fad {

typedef Sacado::Fad::DFad<double> FADDouble;

/**
 * @brief Enumeration <b>Opcode</b> lists the operation codes which can be
 * executed using the execution virtual machine.
 *
 */
enum Opcode {
  OPCODE_UNDEFINED,
  OPCODE_CONSTANT,
  OPCODE_RVALUE,
  OPCODE_STATEMENT,
  OPCODE_ARGUMENT,

  OPCODE_TIERNARY,
  OPCODE_MULTIPLY,
  OPCODE_DIVIDE,
//  OPCODE_MODULUS,
  OPCODE_ADD,
  OPCODE_SUBTRACT,
  OPCODE_UNARY_MINUS,
  OPCODE_FUNCTION,

  OPCODE_EQUAL,
  OPCODE_NOT_EQUAL,
  OPCODE_LESS,
  OPCODE_GREATER,
  OPCODE_LESS_EQUAL,
  OPCODE_GREATER_EQUAL,

  OPCODE_UNARY_NOT,
  OPCODE_LOGICAL_AND,
  OPCODE_LOGICAL_OR,

  OPCODE_EXPONENIATION,
  OPCODE_ASSIGN
};


class Node
{
public:
  explicit Node(Opcode opcode)
    : m_opcode(opcode),
      m_left(0),
      m_right(0),
      m_other(0)
  {}

private:
  explicit Node(const Node &);
  Node &operator=(const Node &);

public:
  ~Node()
  {}

  FADDouble eval() const;

  const Opcode	m_opcode;

  union _data
  {
    struct _constant
    {
      double value;
    } constant;

    struct _variable
    {
      Eval::Variable *variable;
    } variable;

    struct _function
    {
      CFunctionBase *	function;
    } function;
  } m_data;

  Node *		m_left;
  Node *		m_right;
  Node *		m_other;
};


FADDouble
Node::eval() const
{
  switch (m_opcode) {
  case OPCODE_STATEMENT:
    {
      FADDouble value = 0.0;
      for (const Node *statement = this; statement; statement = statement->m_right)
	value = statement->m_left->eval();
      return value;
    }

  case OPCODE_CONSTANT:
    return m_data.constant.value;

  case OPCODE_RVALUE:
    /* Directly access the variable */
    if (m_left)
      return (*m_data.variable.variable)[m_left->eval()];
    else
      return m_data.variable.variable->getValue();

  case OPCODE_MULTIPLY:
    return m_left->eval()*m_right->eval();

  case OPCODE_EXPONENIATION:
    return std::pow(m_left->eval(),m_right->eval());

  case OPCODE_DIVIDE:
    return m_left->eval()/m_right->eval();

//   case OPCODE_MODULUS:
//     return std::fmod(m_left->eval(), m_right->eval());

  case OPCODE_ADD:
    return m_left->eval() + m_right->eval();

  case OPCODE_SUBTRACT:
    return m_left->eval() - m_right->eval();

  case OPCODE_EQUAL:
    return m_left->eval() == m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_NOT_EQUAL:
    return m_left->eval() != m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_LESS:
    return m_left->eval() < m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_GREATER:
    return m_left->eval() > m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_LESS_EQUAL:
    return m_left->eval() <= m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_GREATER_EQUAL:
    return m_left->eval() >= m_right->eval() ? Eval::s_true : Eval::s_false;

  case OPCODE_LOGICAL_AND: {
    FADDouble left = m_left->eval();
    FADDouble right = m_right->eval();
    return  (left != Eval::s_false) && (right != Eval::s_false) ? Eval::s_true : Eval::s_false;
  }
  case OPCODE_LOGICAL_OR: {
    FADDouble left = m_left->eval();
    FADDouble right = m_right->eval();
    
    return (left != Eval::s_false) || (right != Eval::s_false) ? Eval::s_true : Eval::s_false;
  }
  case OPCODE_TIERNARY:
    return m_left->eval() != Eval::s_false ? m_right->eval() : m_other->eval();

  case OPCODE_UNARY_MINUS:
    return -m_right->eval();

  case OPCODE_UNARY_NOT:
    return m_right->eval() == Eval::s_false ? Eval::s_true : Eval::s_false;

  case OPCODE_ASSIGN:
    if (m_left)
      return (*m_data.variable.variable)[m_left->eval()] = m_right->eval();
    else {
      *m_data.variable.variable = m_right->eval();
      return m_data.variable.variable->getValue();
    }

  case OPCODE_FUNCTION:
    {
      FADDouble argv[20];

      int argc = 0;
      for (Node *arg = m_right; arg; arg = arg->m_right)
	argv[argc++] = arg->m_left->eval();

      return (*m_data.function.function)(argc, argv);
    }

  default: // Unknown opcode
    throw std::runtime_error("Evaluation error");
  }
}

namespace Parser {

Node *parseStatements(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseStatement(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseExpression(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseAssign(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator assign, LexemVector::const_iterator to);
Node *parseTerm(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator term, LexemVector::const_iterator to);
Node *parseFactor(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseRelation(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseLogical(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator factor, LexemVector::const_iterator to);
Node *parseUnary(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator unary, LexemVector::const_iterator to);
Node *parseTiernary(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator question, LexemVector::const_iterator colon, LexemVector::const_iterator to);
Node *parseFunction(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lparen, LexemVector::const_iterator rparen, LexemVector::const_iterator to);
Node *parseFunctionArg(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseRValue(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator to);
Node *parseIndex(Eval &eval, LexemVector::const_iterator from, LexemVector::const_iterator lbrack, LexemVector::const_iterator rbrack, LexemVector::const_iterator to);

// Parser productions

Node *
parseStatements(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() == TOKEN_END)
    return NULL;

  if ((*from).getToken() == TOKEN_SEMI)
    return parseStatements(eval, from + 1, to);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  LexemVector::const_iterator it;
  for (it = from; (*it).getToken() != TOKEN_SEMI && (*it).getToken() != TOKEN_END; ++it)
    ;

  Node *statement = eval.newNode(OPCODE_STATEMENT);
  statement->m_left = parseStatement(eval, from, it);

  // Technically, there should be no check for TOKEN_END, but we allow a missing final TOKEN_SEMI
  if ((*it).getToken() != TOKEN_END)
    statement->m_right = parseStatements(eval, it + 1, to);

  return statement;
}


Node *
parseStatement(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  return parseExpression(eval, from, to);
}


Node *
parseExpression(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  int paren_level = 0;					// Paren level
  int brack_level = 0;					// Brack level
  LexemVector::const_iterator lparen_open_it = to;	// First open paren
  LexemVector::const_iterator lparen_close_it = to;	// Corresponding close paren
  LexemVector::const_iterator lbrack_open_it = to;	// First open bracket
  LexemVector::const_iterator lbrack_close_it = to;	// Corresponding close brack
  LexemVector::const_iterator assign_it = to;		// First = at paren_level 0 for assignment
  LexemVector::const_iterator term_it = to;		// Last + or - at paren_level 0 for adding or subtracting
  LexemVector::const_iterator factor_it = to;		// Last * or / at paren_level 0 for multiplying or dividing
  LexemVector::const_iterator relation_it = to;		// Last relational at paren_level 0 for relational operator
  LexemVector::const_iterator logical_it = to;		// Last logical at paren_level 0 for logical operator
  LexemVector::const_iterator question_it = to;		// Last tiernary at paren_level 0 for tiernary operator
  LexemVector::const_iterator colon_it = to;
  LexemVector::const_iterator unary_it = to;		// First +,- at plevel 0 for positive,negative
  LexemVector::const_iterator last_unary_it = to;	// Last +,- found at plevel for for positive,negative

  // Scan the expression for the instances of the above tokens
  for (LexemVector::const_iterator it = from; it != to; ++it) {
    switch((*it).getToken()) {
    case TOKEN_LPAREN:
      if (paren_level == 0 && lparen_open_it == to
	  && brack_level == 0 && lbrack_open_it == to)
	lparen_open_it = it;
      paren_level++;
      break;

    case TOKEN_RPAREN:
      paren_level--;

      if (paren_level == 0 && lparen_close_it == to
	  && brack_level == 0 && lbrack_close_it == to)
	lparen_close_it = it;

      if (paren_level < 0)
	throw std::runtime_error("mismatched parenthesis");
      break;

    case TOKEN_LBRACK:
      if (paren_level == 0 && lparen_open_it == to
	  && brack_level == 0 && lbrack_open_it == to)
	lbrack_open_it = it;
      brack_level++;
      break;

    case TOKEN_RBRACK:
      brack_level--;

      if (paren_level == 0 && lparen_close_it == to
	  && brack_level == 0 && lbrack_close_it == to)
	lbrack_close_it = it;

      if (brack_level < 0)
	throw std::runtime_error("mismatched bracket");
      break;

    case TOKEN_ASSIGN:
      if (paren_level == 0 && assign_it == to)
	assign_it = it;
      break;

    case TOKEN_QUESTION:
      if (paren_level == 0 && question_it == to)
	question_it = it;
      break;

    case TOKEN_COLON:
      if (paren_level == 0) // && colon_it == to)
	colon_it = it;
      break;

    case TOKEN_EXPONENTIATION:
    case TOKEN_MULTIPLY:
    case TOKEN_DIVIDE:
//    case TOKEN_PERCENT:
      if (paren_level == 0) // && factor_it == to)
	factor_it = it;
      break;

    case TOKEN_EQUAL:
    case TOKEN_NOT_EQUAL:
    case TOKEN_LESS:
    case TOKEN_GREATER:
    case TOKEN_LESS_EQUAL:
    case TOKEN_GREATER_EQUAL:
      if (paren_level == 0 && relation_it == to)
	relation_it = it;
      break;

    case TOKEN_LOGICAL_AND:
    case TOKEN_LOGICAL_OR:
      if (paren_level == 0 && logical_it == to)
	logical_it = it;
      break;

    case TOKEN_PLUS:
    case TOKEN_MINUS:
      if (paren_level == 0) {
	// After any of these, we are a unary operator, not a term
	if (it == from || it == assign_it + 1
	    || it == term_it + 1 || it == factor_it + 1
	    || it == last_unary_it + 1)
	{ // Unary operator
	  if (unary_it == to) // First unary operator?
	    unary_it = it;
	  last_unary_it = it;
	}
	else { // Term
	  term_it = it;
	}
      }
      break;

    case TOKEN_NOT:
      if (paren_level == 0) {
	if (unary_it == to) /// First unary operator
	  unary_it = it;
	last_unary_it = it;
      }
      break;

    default:
      break;
    }
  }

  if (paren_level != 0) // paren_level should now be zero */
    throw std::runtime_error("mismatched parenthesis");

  // This implement the operator hiearchy
  // Assignment
  if (assign_it != to)
    return parseAssign(eval, from, assign_it, to);

  // Tiernary operator
  if (question_it != to || colon_it != to)
    return parseTiernary(eval, from, question_it, colon_it, to);

  // Logical
  if (logical_it != to)
    return parseLogical(eval, from, logical_it, to);

  // Relational
  if (relation_it != to)
    return parseRelation(eval, from, relation_it, to);

  // Term
  if (term_it != to)
    return parseTerm(eval, from, term_it, to);

  // Factor
  if (factor_it != to)
    return parseFactor(eval, from, factor_it, to);

  // Unary
  if (unary_it != to)
    return parseUnary(eval, from, unary_it, to);

  // Parenthetical
  if (lparen_open_it != to) {
    if (lparen_open_it == from) {
      if (lparen_close_it == to - 1 && lparen_close_it - lparen_open_it > 1)
	return parseExpression(eval, lparen_open_it + 1, lparen_close_it);
      else
	throw std::runtime_error("syntax error parsing parentheses");
    }

    // Function
    if (lparen_open_it == from + 1) {
      if (lparen_close_it == to - 1)
	return parseFunction(eval, from, lparen_open_it, lparen_close_it, to);
      else // Closing paren not at to
	throw std::runtime_error("syntax error 2");
    }

    throw std::runtime_error("syntax error 3");
  }

  // Bracket
  if (lbrack_open_it != to) {

    // Array index
    if (lbrack_open_it == from + 1) {
      if (lbrack_close_it == to - 1)
	return parseIndex(eval, from, lbrack_open_it, lbrack_close_it, to);
      else // Closing brack not at to
	throw std::runtime_error("syntax error 2");
    }

    throw std::runtime_error("syntax error 3");
  }

  // R-Value
  return parseRValue(eval, from, to);
}


Node *
parseAssign(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	assign_it,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER) //  || from + 1 != assign_it) {
    throw std::runtime_error("stk_classic::expreval::parseAssign: expected identifier");

  Node *assign;

  if ((*(from + 1)).getToken() == TOKEN_ASSIGN || (*(from + 1)).getToken() == TOKEN_LBRACK) {
    assign = eval.newNode(OPCODE_ASSIGN);

    assign->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
    assign->m_right = parseExpression(eval, assign_it + 1, to);

    if ((*(from + 1)).getToken() == TOKEN_LBRACK)
      assign->m_left = parseExpression(eval, from + 2, assign_it - 1);
  }
  else
    throw std::runtime_error("syntax error");

  return assign;
}


Node *
parseTerm(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	term_it,
  LexemVector::const_iterator	to)
{
  Node *term = eval.newNode((*term_it).getToken() == TOKEN_PLUS ? OPCODE_ADD : OPCODE_SUBTRACT);

  term->m_left = parseExpression(eval, from, term_it);
  term->m_right = parseExpression(eval, term_it + 1, to);

  return term;
}


Node *
parseFactor(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	factor_it,
  LexemVector::const_iterator	to)
{
  Node *factor = eval.newNode((*factor_it).getToken() == TOKEN_MULTIPLY ? OPCODE_MULTIPLY : OPCODE_DIVIDE);

  factor->m_left = parseExpression(eval, from, factor_it);
  factor->m_right = parseExpression(eval, factor_it + 1, to);

  return factor;
}


Node *
parseRelation(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	relation_it,
  LexemVector::const_iterator	to)
{
  Opcode relation_opcode = OPCODE_UNDEFINED;

  switch ((*relation_it).getToken()) {
  case TOKEN_EQUAL:
    relation_opcode = OPCODE_EQUAL;
    break;

  case TOKEN_NOT_EQUAL:
    relation_opcode = OPCODE_NOT_EQUAL;
    break;

  case TOKEN_LESS:
    relation_opcode = OPCODE_LESS;
    break;

  case TOKEN_GREATER:
    relation_opcode = OPCODE_GREATER;
    break;

  case TOKEN_LESS_EQUAL:
    relation_opcode = OPCODE_LESS_EQUAL;
    break;

  case TOKEN_GREATER_EQUAL:
    relation_opcode = OPCODE_GREATER_EQUAL;
    break;

  default:
    break;
  }

  Node *relation = eval.newNode(relation_opcode);

  relation->m_left = parseExpression(eval, from, relation_it);
  relation->m_right = parseExpression(eval, relation_it + 1, to);

  return relation;
}


Node *
parseLogical(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	logical_it,
  LexemVector::const_iterator	to)
{
  Node *logical = eval.newNode(((*logical_it).getToken() == TOKEN_LOGICAL_AND ? OPCODE_LOGICAL_AND : OPCODE_LOGICAL_OR));

  logical->m_left = parseExpression(eval, from, logical_it);
  logical->m_right = parseExpression(eval, logical_it + 1, to);

  return logical;
}


Node *
parseTiernary(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	question_it,
  LexemVector::const_iterator	colon_it,
  LexemVector::const_iterator	to)
{
  if (question_it == to || colon_it == to)
    throw std::runtime_error("syntax error parsing ?: operator");

  Node *tiernary = eval.newNode(OPCODE_TIERNARY);

  tiernary->m_left = parseExpression(eval, from, question_it);
  tiernary->m_right = parseExpression(eval, question_it + 1, colon_it);
  tiernary->m_other = parseExpression(eval, colon_it + 1, to);

  return tiernary;
}


Node *
parseUnary(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	unary_it,
  LexemVector::const_iterator	to)
{
  /* If it is a positive, just parse the internal of it */
  if ((*unary_it).getToken() == TOKEN_PLUS)
    return parseExpression(eval, unary_it + 1, to);
  else if ((*unary_it).getToken() == TOKEN_MINUS) {
    Node *unary = eval.newNode(OPCODE_UNARY_MINUS);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else if ((*unary_it).getToken() == TOKEN_NOT) {
    Node *unary = eval.newNode(OPCODE_UNARY_NOT);
    unary->m_right = parseExpression(eval, unary_it + 1, to);
    return unary;
  }
  else
    throw std::runtime_error("syntax error parsing unary operator");
}


Node *
parseFunction(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	lparen,
  LexemVector::const_iterator	rparen,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing function");

  const std::string &function_name = (*from).getString();

  CFunctionBase *c_function = NULL;
  Eval::CFunctionMap::iterator it = Eval::getCFunctionMap().find(function_name);
  if (it != Eval::getCFunctionMap().end())
    c_function = (*it).second;

//   if (!c_function)
//     throw std::runtime_error(std::string("Undefined function ") + function_name);

  Node *function = eval.newNode(OPCODE_FUNCTION);
  function->m_data.function.function = c_function;

  if (!c_function)
    eval.getUndefinedFunctionSet().insert(function_name);

  function->m_right = parseFunctionArg(eval, lparen + 1, rparen);

  //   if (!c_function)
  //     throw std::runtime_error(std::string("Undefined function ") + function_name);

  return function;
}


Node *
parseIndex(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	lbrack,
  LexemVector::const_iterator	rbrack,
  LexemVector::const_iterator	to)
{
  if ((*from).getToken() != TOKEN_IDENTIFIER)
    throw std::runtime_error("syntax error parsing array");

  Node *index = eval.newNode(OPCODE_RVALUE);
  index->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
  index->m_left = parseExpression(eval, lbrack + 1, rbrack);

  return index;
}


Node *
parseFunctionArg(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if (from == to)
    return NULL;

  LexemVector::const_iterator it;
  for (it = from; it != to && (*it).getToken() != TOKEN_COMMA; ++it)
    ;

  Node *argument = eval.newNode(OPCODE_ARGUMENT);
  argument->m_left = parseExpression(eval, from, it);
  if (it != to)
    argument->m_right = parseFunctionArg(eval, it + 1, to);
  return argument;
}


Node *
parseRValue(
  Eval &			eval,
  LexemVector::const_iterator	from,
  LexemVector::const_iterator	to)
{
  if (from + 1 != to)
    throw std::runtime_error(std::string("r-value not allowed following ") + (*from).getString());

  switch ((*from).getToken()) {
  case TOKEN_IDENTIFIER:
    {
      Eval::ConstantMap::iterator it = Eval::getConstantMap().find((*from).getString());
      if (it != Eval::getConstantMap().end()) {
	Node *constant = eval.newNode(OPCODE_CONSTANT);
	constant->m_data.constant.value = (*it).second;
	return constant;
      }

      // Define a variable
      else {
	Node *variable = eval.newNode(OPCODE_RVALUE);
	variable->m_data.variable.variable = eval.getVariableMap()[(*from).getString()];
	return variable;
      }
    }

  case TOKEN_REAL_CONSTANT:
  case TOKEN_INTEGER_CONSTANT:
    {
      Node *constant = eval.newNode(OPCODE_CONSTANT);
      constant->m_data.constant.value = (*from).getValue<double>();
      return constant;
    }

  default:
    throw std::runtime_error("invalid rvalue");
  }
}

} // namespace Parser

extern "C" {
  typedef FADDouble (*CExtern0)();
  typedef FADDouble (*CExtern1)(FADDouble);
  typedef FADDouble (*CExtern2)(FADDouble, FADDouble);
}


template <>
class CFunction<CExtern0> : public CFunctionBase
{
public:
  typedef CExtern0 Signature;

  explicit CFunction<CExtern0>(Signature function)
    : CFunctionBase(0),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)();
  }

private:
  Signature	m_function;
};


template <>
class CFunction<CExtern1> : public CFunctionBase
{
public:
  typedef CExtern1 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(1),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)(argv[0]);
  }

private:
  Signature	m_function;
};


template <>
class CFunction<CExtern2> : public CFunctionBase
{
public:
  typedef CExtern2 Signature;

  explicit CFunction<Signature>(Signature function)
    : CFunctionBase(2),
      m_function(function)
  {}

  virtual ~CFunction()
  {}

  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    if (argc != getArgCount())
      throw std::runtime_error("Argument count mismatch");

    return (*m_function)(argv[0], argv[1]);
  }

private:
  Signature	m_function;
};


typedef CFunction<CExtern0> CFunction0;
typedef CFunction<CExtern1> CFunction1;
typedef CFunction<CExtern2> CFunction2;


extern "C" {
  static FADDouble real_rand()  {
    return (FADDouble) std::rand() / ((FADDouble)(RAND_MAX) + 1.0);
  }

  static FADDouble real_srand(FADDouble x)  {
    std::srand(0); // (int) x);
    return 0.0;
  }

  static FADDouble randomize()  {
#ifndef SIERRA_SRAND_OK
    std::srand((::clock() + 1024) * ::time(NULL));
#endif
    return 0.0;
  }
}


const double
Eval::s_e	= 2.7182818284590452354;

const double
Eval::s_pi	= 3.14159265358979323846;

const double
Eval::s_false	= 0.0;

const double
Eval::s_true	= 1.0;


Eval::VariableMap::Resolver &
Eval::VariableMap::getDefaultResolver()
{
  static DefaultResolver default_resolver;

  return default_resolver;
}


Eval::ConstantMap &
Eval::getConstantMap()
{
  static Eval::ConstantMap s_constantMap;

  if (s_constantMap.empty()) {
    s_constantMap["E"] = s_e;
    s_constantMap["PI"] = s_pi;
    s_constantMap["FALSE"] = Eval::s_false;
    s_constantMap["TRUE"] = Eval::s_true;
  }

  return s_constantMap;
}


class CFunctionExp : public CFunctionBase
{
public:
  CFunctionExp()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return exp(argv[0]);
  }
};

  
class CFunctionLog : public CFunctionBase
{
public:
  CFunctionLog()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return log(argv[0]);
  }
};

  
class CFunctionLog10 : public CFunctionBase
{
public:
  CFunctionLog10()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return log10(argv[0]);
  }
};

  
class CFunctionPow : public CFunctionBase
{
public:
  CFunctionPow()
    : CFunctionBase(2)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return pow(argv[0], argv[1]);
  }
};

  
class CFunctionSqrt : public CFunctionBase
{
public:
  CFunctionSqrt()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return sqrt(argv[0]);
  }
};

  
class CFunctionAcos : public CFunctionBase
{
public:
  CFunctionAcos()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return acos(argv[0]);
  }
};

  
class CFunctionAsin : public CFunctionBase
{
public:
  CFunctionAsin()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return asin(argv[0]);
  }
};

  
class CFunctionAtan : public CFunctionBase
{
public:
  CFunctionAtan()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return atan(argv[0]);
  }
};

  
// class CFunctionAtan2 : public CFunctionBase
// {
// public:
//   CFunctionAtan2()
//     : CFunctionBase(2)
//   {}
  
//   virtual FADDouble operator()(int argc, const FADDouble *argv) {
//     return atan2(argv[0], argv[1]);
//   }
// };

  
class CFunctionCos : public CFunctionBase
{
public:
  CFunctionCos()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return cos(argv[0]);
  }
};

  
class CFunctionCosh : public CFunctionBase
{
public:
  CFunctionCosh()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return cosh(argv[0]);
  }
};

  
class CFunctionSin : public CFunctionBase
{
public:
  CFunctionSin()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return sin(argv[0]);
  }
};

  
class CFunctionSinh : public CFunctionBase
{
public:
  CFunctionSinh()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return sinh(argv[0]);
  }
};

  
class CFunctionTan : public CFunctionBase
{
public:
  CFunctionTan()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return tan(argv[0]);
  }
};

  
class CFunctionTanh : public CFunctionBase
{
public:
  CFunctionTanh()
    : CFunctionBase(1)
  {}
  
  virtual FADDouble operator()(int argc, const FADDouble *argv) {
    return tanh(argv[0]);
  }
};

  
Eval::CFunctionMap::CFunctionMap() 
{
  (*this)["rand"] = new CFunction0(real_rand);
  (*this)["srand"] = new CFunction1(real_srand);
  (*this)["randomize"] = new CFunction0(randomize);

  (*this)["exp"] = new CFunctionExp();
  (*this)["ln"] = new CFunctionLog();
  (*this)["log"] = new CFunctionLog();
  (*this)["log10"] = new CFunctionLog10();
  (*this)["pow"] = new CFunctionPow();
  (*this)["sqrt"] = new CFunctionSqrt();

  (*this)["acos"] = new CFunctionAcos();
  (*this)["asin"] = new CFunctionAsin();
  (*this)["atan"] = new CFunctionAtan();
  (*this)["cos"] = new CFunctionCos();
  (*this)["cosh"] = new CFunctionCosh();
  (*this)["sin"] = new CFunctionSin();
  (*this)["sinh"] = new CFunctionSinh();
  (*this)["tan"] = new CFunctionTan();
  (*this)["tanh"] = new CFunctionTanh();
}

Eval::CFunctionMap::~CFunctionMap()
{
  for (CFunctionMap::iterator it = begin(); it != end(); ++it)
    delete (*it).second;
}


Eval::CFunctionMap &
Eval::getCFunctionMap()
{
  static CFunctionMap s_functionMap;

  return s_functionMap;
}

Eval::Eval(
  VariableMap::Resolver &		resolver,
  const std::string &			expression)
  : m_variableMap(resolver),
    m_expression(expression),
    m_parseStatus(false),
    m_headNode(0)
{}


Eval::~Eval()
{
  for (std::vector<Node *>::iterator it = m_nodes.begin(); it != m_nodes.end(); ++it)
    delete (*it);
}


Node *
Eval::newNode(
  int           opcode)
{
  Node *new_node = new Node((Opcode) opcode);
  m_nodes.push_back(new_node);
  return new_node;
}


void
Eval::syntax()
{
  m_syntaxStatus = false;
  m_parseStatus = false;

  try {
    /* Validate the characters */
    LexemVector lex_vector = tokenize(m_expression);

    /* Call the multiparse routine to parse subexpressions */
    m_headNode = Parser::parseStatements(*this, lex_vector.begin(), lex_vector.end());

    m_syntaxStatus = true;
  }
  catch (std::runtime_error & /* x */) {
//     x << " while parsing expression: " << m_expression;
//     RuntimeDoomed() << x.what();
    throw;
  }
}


void
Eval::parse()
{
  try {
    syntax();

    if (m_syntaxStatus) {
      if (!m_undefinedFunctionSet.empty()) {
        std::ostringstream strout;
        strout << "In expression '" << m_expression << "', the following functions are not defined:" << std::endl;
        //	for (iteration<UndefinedFunctionSet> it(m_undefinedFunctionSet); it; ++it)
        for (UndefinedFunctionSet::iterator it = m_undefinedFunctionSet.begin(); it != m_undefinedFunctionSet.end(); ++it)
          strout << (*it) << std::endl;
        throw std::runtime_error(strout.str());
      }

      resolve();

      m_parseStatus = true;
    }
  }
  catch (std::runtime_error & /* x */) {
//     x << " while parsing expression: " << m_expression;
//     RuntimeDoomed() << x.what();
    throw;
  }
}


void
Eval::resolve()
{
  for (VariableMap::iterator it = m_variableMap.begin(); it != m_variableMap.end(); ++it) {
    m_variableMap.getResolver().resolve(it);
  }
}


FADDouble
Eval::evaluate() const
{
  /* Make sure it was parsed successfully */
  if (!m_parseStatus)
    throw std::runtime_error(std::string("Expression '") + m_expression + "' did not parse successfully");

  return m_headNode->eval();
}

} // namespace fad
} // namespace expreval
} // namespace stk_classic
