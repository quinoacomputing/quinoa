//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon Oct  7 11:41:06 2013
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
#include <Quinoa/CmdLine/Keywords.h>
#include <Quinoa/InputDeck/InputDeck.h>

namespace quinoa {
//! Grammar definition: state, actions, grammar
namespace cmd {

  using namespace pegtl;
  using namespace tk::grm;

  // State

  //! Everything is stored in Stack during parsing
  using Stack = InputDeck;

  // Actions

  //! convert and put value in state at position given by tags
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store<tags...>(value);
    }
  };

  // Grammar

  //! match verbose cmdline keyword
  template< class keyword >
  struct verbose :
         seq< string<'-','-'>, typename keyword::pegtl_string, space > {};

  //! match alias cmdline keyword
  template< class keyword >
  struct alias :
         seq< one<'-'>,
              typename keyword::pegtl_alias,
              sor<space, unknown<Stack,Error::ALIAS> > > {};

  //! read 'keyword' in either verbose or alias form
  template< class keyword >
  struct readkw :
         sor< verbose<keyword>, alias<keyword> > {};

  //! parse input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions... >, blank, space > {};

  //! process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = any >
  struct process :
         ifmust< readkw<keyword>,
                 parse< sor<keywords, apply<error<Stack,Error::MISSING>>>,
                        insert> > {};

  //! control (i.e., input deck) file
  struct control :
         process<kw::control, store<ctr::io,ctr::control>> {};

  //! input file
  struct input :
         process<kw::input, store<ctr::io,ctr::input>> {};

  //! output file
  struct output :
         process<kw::output, store<ctr::io,ctr::output>> {};

  //! pdf output file
  struct pdf :
         process<kw::pdf, store<ctr::io,ctr::pdf>> {};

  //! glob output file
  struct glob :
         process<kw::glob, store<ctr::io,ctr::glob>> {};

  //! stat output file
  struct stat :
         process<kw::stat, store<ctr::io,ctr::stat>> {};

  //! command line keywords
  struct keywords :
         sor< control,
              input,
              output,
              pdf,
              glob,
              stat > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         until< eof, sor<keywords, unknown<Stack,Error::KEYWORD>> > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
