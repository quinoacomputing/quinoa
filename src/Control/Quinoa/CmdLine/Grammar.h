//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Thu Oct  3 12:19:48 2013
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
#include <Quinoa/CmdLine/Keywords.h>

namespace quinoa {
//! Grammar definition: state, actions, grammar
namespace cmd {

  using namespace pegtl;

  // State

  //! Everything is stored in Stack during parsing
  using Stack = InputDeck;

  //! Command line parser error types
  enum class Error : uint8_t { KEYWORD,
                               ALIAS,
                               MISSING };

  static const std::map<Error, std::string> err_msg( {
    { Error::KEYWORD, "Unknown keyword" },
    { Error::ALIAS, "Alias keyword too long" },
    { Error::MISSING, "Required filename missing" }
  } );

  // Actions

  //! error handler
  template< Error key >
  struct error : action_base< error<key> > {
    static void apply(const std::string& value, Stack& stack) {
      const auto& msg = err_msg.find(key);
      if (msg != err_msg.end()) {
        if (!value.empty()) {
          Throw(ExceptType::FATAL,
                "Error while parsing '" + value + "' in the command line. " +
                msg->second + ".");
        } else {
          Throw(ExceptType::FATAL,
                "Error while parsing the command line. " + msg->second + ".");
        }
      } else {
        Throw(ExceptType::FATAL, "Unknown command line parser error.");
      }
      IGNORE(stack);
    }
  };

  //! convert and put value in state at position given by tags
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store<tags...>(value);
    }
  };

  // Grammar

  //! read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  // match unknown keyword and handle error
  template< Error key >
  struct unknown :
         pad< ifapply< trim<any, space>, error<key> >, blank, space > {};

  //! match verbose cmdline keyword
  template< class keyword >
  struct verbose :
         seq< string<'-','-'>, typename keyword::pegtl_string, space > {};

  //! match alias cmdline keyword
  template< class keyword >
  struct alias :
         seq< one<'-'>,
              typename keyword::pegtl_alias,
              sor<space, unknown<Error::ALIAS> > > {};

  //! read 'keyword' in either verbose or alias form
  template< class keyword >
  struct read :
         sor< verbose<keyword>, alias<keyword> > {};

  //! parse input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions... >, blank, space > {};

  //! process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = any >
  struct process :
         ifmust< read<keyword>,
                 parse< sor<keywords, apply<error<Error::MISSING>>>,
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
         until< eof, sor<keywords, unknown<Error::KEYWORD>> > {};

} // cmd::
} // quinoa::

#endif // QuinoaCmdLineGrammar_h
