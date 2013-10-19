//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Fri Oct 18 12:53:17 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

  using namespace pegtl;
  using namespace tk::grm;

  //! PEGTLParsed type specialized to RNGTest's command line parser
  using PEGTLCmdLine =
    quinoa::ctr::PEGTLParsed< ctr::CmdLine, string_input<ctr::Location> >;

  // RNGTest's CmdLine state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLCmdLine;

  // RNGTest's CmdLine actions

  // RNGTest's CmdLine grammar

  //! control (i.e., input deck) file
  struct control :
         process_cmd<Stack, kw::control, Store<Stack,ctr::io,ctr::control>> {};

  //! command line keywords
  struct keywords :
         sor< control > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         until< eof, sor<keywords, unknown<Stack,Error::KEYWORD>> > {};

} // cmd::
} // rngtest::

#endif // RNGTestCmdLineGrammar_h
