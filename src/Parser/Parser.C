//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 08:26:35 PM MST
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

  // Keywords

  #include <Keywords.def.h>

  // State

  // Actions

  struct do_comment : action_base< do_comment > {
    static void apply(const std::string& m) {
      std::cout << "COMMENT: \"" << m << "\"" << endl;
    }
  };

  struct do_keyword : action_base< do_keyword > {
    static void apply(const std::string& m) {
      std::cout << "KEYWORD: \"" << m << "\"" << endl;
    }
  };

  struct do_input : action_base< do_input > {
    static void apply(const std::string& m) {
      std::cout << "INPUT  : \"" << m << "\"" << endl;
    }
  };

  // Grammar

  struct trim_input :
         seq< alnum, until< at<space> > > {};

  struct input :
         seq< star<blank>, ifapply< trim_input, do_input>, space > {};

  struct trim_keyw :
         seq< keyword::any, until< at<space> > > {};

  struct keyw :
         seq< star<blank>, ifapply< trim_keyw, do_keyword>, space > {};

  struct token :
         sor< keyw, input > {};

  struct trim_comment :
         seq< one<'#'>, until< at<eol> > > {};

  struct comment :
         seq< star<blank>, ifapply< trim_comment, do_comment >, eol> {};
         //seq< one<'#'>, until<at<eol>>, eol> {};

  struct read_file :
         until< eof, sor<token, comment> > {};

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

  cout << "==== PARSE START ====" << endl;
  pegtl::basic_parse_file< grammar::read_file >( m_filename );
  cout << "==== PARSE END ====" << endl << endl;
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
