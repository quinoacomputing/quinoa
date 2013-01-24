//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Wed 23 Jan 2013 11:02:39 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Parser base
  \details   Parser base
*/
//******************************************************************************

#include <pegtl.hh>

#include <Parser.h>
#include <IOException.h>

namespace grammar {

  using namespace pegtl;
  using namespace pegtl::ascii;

  // State

  typedef int value_type;
  typedef vector< value_type > stack_type;

  // Keywords

  //typedef pegtl::string<t,i,t,l,e> keyword_title;
  struct keyword_title :
         seq< one<'T','t'>,
              one<'I','i'>,
              one<'T','t'>,
              one<'L','l'>,
              one<'E','e'> > {};
  //typedef pegtl::string<S,P,I,N,S,F,l,o,w> keyword_SPINSFlow;

  // Actions

  struct comment_action : action_base< comment_action > {
    static void apply(const std::string& m) {
      std::cout << "COMMENT: " << m;
    }
  };

  struct rule_action : action_base< rule_action > {
    static void apply(const std::string& m) {
      std::cout << "RULE   : " << m;
    }
  };

  struct title_action : action_base< title_action > {
    static void apply(const std::string& m) {
      std::cout << "TITLE  : " << m;
    }
  };

  // Grammar

  struct read_comment :
         ifapply< seq< until<one<'#'>>, until<eol> >, comment_action > {};

  struct read_rule :
	 ifapply< seq< until<alpha>, until<eol> >, rule_action > {};

  struct read_title :
	 ifapply< seq< until<keyword_title>, until<eol> >, title_action > {};

  struct read_line :
         sor< read_comment, read_title, read_rule > {};

  struct read_file :
         until< eof, read_line > {};

} // namespace grammar

using namespace Quinoa;

Parser::Parser(const string& filename) : m_filename(filename)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
//  m_q.open(m_filename, ifstream::in);
//  Assert(m_q.good(), IOException,FATAL,IO_FAILED_OPEN,m_filename);

  cout << "==== PARSING START ====" << endl;
  pegtl::basic_parse_file< grammar::read_file >( m_filename );
  cout << "==== PARSING END ====" << endl << endl;
}

Parser::~Parser()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
//  m_q.close();

  // No exception leaves a destructor: if the above close() fails, we only emit
  // a warning, thus we avoid terminate if an exception is propagating through.
//  if (m_q.fail())
//    cout << "WARNING: Failed to close file: " << m_filename << endl;
}
