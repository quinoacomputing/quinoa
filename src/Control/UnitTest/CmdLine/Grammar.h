//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:21:59 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
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
#include <UnitTest/CmdLine/Keywords.h>

namespace unittest {
//! Grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to UnitTest's command line parser
  using PEGTLCmdLine =
    quinoa::ctr::PEGTLParsed< ctr::CmdLine,
                              pegtl::string_input< ctr::Location > >;

  // UnitTest's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // UnitTest's CmdLine actions

  // UnitTest's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      tk::kw::verbose,
                                      tk::tag::verbose > {};

  //! io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< Stack,
                               keyword,
                               tk::grm::Store< Stack, tag::io, io_tag > > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     io< kw::input, tag::input >,
                     io< kw::output, tag::output > > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // unittest::

#endif // UnitTestCmdLineGrammar_h
