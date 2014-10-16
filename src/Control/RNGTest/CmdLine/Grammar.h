//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 10:10:16 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     RNGTest's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef RNGTestCmdLineGrammar_h
#define RNGTestCmdLineGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <Grammar.h>
#include <PEGTLParsed.h>
#include <RNGTest/CmdLine/Keywords.h>

namespace rngtest {
//! Grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to RNGTest's command line parser
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  // RNGTest's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // RNGTest's CmdLine actions

  // RNGTest's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      tk::kw::verbose,
                                      tk::tag::verbose > {};

  //! control (i.e., input deck) file name
  struct control :
         tk::grm::process_cmd< Stack,
                               kw::control,
                               tk::grm::Store< Stack,
                                               tag::io,
                                               tag::control > > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< verbose, control > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // rngtest::

#endif // RNGTestCmdLineGrammar_h
