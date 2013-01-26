//******************************************************************************
/*!
  \file      src/Parser/Parser.C
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 11:15:15 PM MST
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

  struct do_input : action_base< do_input > {
    static void apply(const std::string& m) {
      std::cout << "INPUT  : \"" << m << "\"" << endl;
    }
  };

  struct parse_title : action_base< parse_title > {
    static void apply(const std::string& m) {
      std::cout << "TITLE  : \"" << m << "\"" << endl;
    }
  };

  // Utilities

  template< class what, class erased >
  struct trim :
         seq< what, until< at<erased> > > {};

  template< class what >
  struct read :
         pad< trim<what, space>, blank, space > {};

  // Grammar

  struct quoted :
         trim< not_one<'"'>, one<'"'> > {};

  struct read_title :
         ifmust< one<'"'>, ifapply< quoted, parse_title >, one<'"'>, space > {};

  struct process_title :
         ifmust< read<keyword::title>, read_title > {};

  struct keyword_main :
         sor< process_title,
              read< keyword::end > > {};

  struct keyword_physics :
         sor< read< keyword::HomogeneousDirichlet >,
              read< keyword::HomogeneousGeneralizedDirichlet >,
              read< keyword::SPINSFlow > > {};

  struct keyword_any :
         sor< keyword_main,
              keyword_physics > {};

  struct input :
         pad< ifapply< trim<alnum, space>, do_input >, blank, space > {};

  struct comment :
         //pad< ifapply< trim<one<'#'>,eol>, do_comment >, blank, eol > {};
         pad< trim<one<'#'>,eol>, blank, eol> {};

  struct ignore :
         sor< comment, eol > {};

  struct read_file :
         until< eof, sor<keyword_any, input, ignore> > {};

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
