//******************************************************************************
/*!
  \file      src/Control/Walker/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:28:33 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef WalkerCmdLineGrammar_h
#define WalkerCmdLineGrammar_h

#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Keywords.h>

namespace walker {
//! Command line parser grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to Walker's command line parser
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  // Walker's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // Walker's CmdLine actions

  // Walker's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      kw::verbose,
                                      tag::verbose > {};

  //! virtualization parameter
  struct virtualization :
         tk::grm::process_cmd< Stack,
                               kw::virtualization,
                               tk::grm::Store< Stack, tag::virtualization >,
                               tk::grm::number > {};

  //! io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< Stack,
                               keyword,
                               tk::grm::Store< Stack, tag::io, io_tag > > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< verbose,
                     virtualization,
                     io< kw::control, tag::control >,
                     io< kw::pdf, tag::pdf >,
                     io< kw::stat, tag::stat > > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // walker::

#endif // WalkerCmdLineGrammar_h
