//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:32:26 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     UnitTest's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef UnitTestCmdLineGrammar_h
#define UnitTestCmdLineGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Keywords.h>

namespace unittest {
//! Grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to UnitTest's command line parser
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  // UnitTest's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // UnitTest's CmdLine actions

  // UnitTest's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      kw::verbose,
                                      tag::verbose > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< verbose > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // unittest::

#endif // UnitTestCmdLineGrammar_h
