//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Grammar.h
  \author    J. Bakosi
  \date      Mon Oct  7 14:32:34 2013
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
#include <RNGTest/CmdLine/Keywords.h>
#include <RNGTest/InputDeck/InputDeck.h>

namespace rngtest {
//! Grammar definition: state, actions, grammar
namespace cmd {

  using namespace pegtl;

  // State

  //! Everything is stored in Stack during parsing
  using Stack = ctr::InputDeck;

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
          Throw(tk::ExceptType::FATAL,
                "Error while parsing '" + value + "' in the command line. " +
                msg->second + ".");
        } else {
          Throw(tk::ExceptType::FATAL,
                "Error while parsing the command line. " + msg->second + ".");
        }
      } else {
        Throw(tk::ExceptType::FATAL, "Unknown command line parser error.");
      }
      IGNORE(stack);    // suppress compiler warning: parameter never referenced
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
                 parse< sor<keywords, apply<error<Error::MISSING>>>,
                        insert> > {};

  //! control (i.e., input deck) file
  struct control :
         process<kw::control, store<ctr::io,ctr::control>> {};

  //! command line keywords
  struct keywords :
         sor< control > {};

  //! entry point: parse keywords and until end of string
  struct read_string :
         until< eof, sor<keywords, unknown<Error::KEYWORD>> > {};

} // cmd::
} // rngtest::

#endif // RNGTestCmdLineGrammar_h
