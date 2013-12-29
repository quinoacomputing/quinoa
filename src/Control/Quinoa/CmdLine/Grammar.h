//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Sat 28 Dec 2013 07:42:07 PM MST
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

  //! io parameter
  template< typename keyword, typename io_tag >
  struct io :
         tk::grm::process_cmd< Stack,
                               keyword,
                               tk::grm::Store< Stack, ctr::io, io_tag > > {};

  //! command line keywords
  struct keywords :
         pegtl::sor< io< kw::control, ctr::control >,
                     io< kw::input, ctr::input >,
                     io< kw::output, ctr::output >,
                     io< kw::pdf, ctr::pdf >,
                     io< kw::glob, ctr::glob >,
                     io< kw::stat, ctr::stat > > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         tk::grm::read_string< Stack, keywords > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
