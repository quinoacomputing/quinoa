//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon 11 Nov 2013 09:32:03 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line grammar definition
  \details   Grammar definition for parsing the command line. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef QuinoaCmdLineGrammar_h
#define QuinoaCmdLineGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Quinoa/CmdLine/Keywords.h>

namespace quinoa {
//! Command line parser grammar definition: state, actions, grammar
namespace cmd {

  //! PEGTLParsed type specialized to Quinoa's command line parser
  using PEGTLCmdLine = ctr::PEGTLParsed< ctr::CmdLine,
                                         pegtl::string_input< ctr::Location > >;

  // Quinoa's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // Quinoa's CmdLine actions

  // Quinoa's CmdLine grammar

  //! control (i.e., input deck) file
  struct control :
         tk::grm::process_cmd< Stack,
                               kw::control,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::control > > {};

  //! input file
  struct input :
         tk::grm::process_cmd< Stack,
                               kw::input,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::input > > {};

  //! output file
  struct output :
         tk::grm::process_cmd< Stack,
                               kw::output,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::output > > {};

  //! pdf output file
  struct pdf :
         tk::grm::process_cmd< Stack,
                               kw::pdf,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::pdf > > {};

  //! glob output file
  struct glob :
         tk::grm::process_cmd< Stack,
                               kw::glob,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::glob > > {};

  //! stat output file
  struct stat :
         tk::grm::process_cmd< Stack,
                               kw::stat,
                               tk::grm::Store< Stack,
                                               ctr::io,
                                               ctr::stat > > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< control,
                     input,
                     output,
                     pdf,
                     glob,
                     stat > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
