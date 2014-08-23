//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Fri 22 Aug 2014 10:47:02 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef QuinoaCmdLineGrammar_h
#define QuinoaCmdLineGrammar_h

#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Quinoa/CmdLine/Keywords.h>

namespace quinoa {
//! Command line parser grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to Quinoa's command line parser
  using PEGTLCmdLine =
    tk::ctr::PEGTLParsed< ctr::CmdLine, pegtl::string_input< ctr::Location > >;

  // Quinoa's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // Quinoa's CmdLine actions

  // Quinoa's CmdLine grammar

  //! verbose (i.e., verbose or quiet output)
  struct verbose :
         tk::grm::process_cmd_switch< Stack,
                                      tk::kw::verbose,
                                      tk::tag::verbose > {};

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
                     io< kw::input, tag::input >,
                     io< kw::output, tag::output >,
                     io< kw::pdf, tag::pdf >,
                     io< kw::glob, tag::glob >,
                     io< kw::stat, tag::stat > > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
