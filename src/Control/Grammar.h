//******************************************************************************
/*!
  \file      src/Control/Grammar.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 11:07:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common of grammars
  \details   Common of grammars
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <Macro.h>
#include <Exception.h>
#include <Option.h>

namespace quinoa {
//! Grammar definition: state, actions, grammar
namespace grm {

  using namespace pegtl;

  // Common auxiliary functions

  //! Parser error types
  enum class Error : uint8_t { KEYWORD,
                               MOMENT,
                               QUOTED,
                               LIST };

  static const std::map<Error, std::string> err_msg( {
    { Error::KEYWORD, "Unknown keyword" },
    { Error::MOMENT, "Unknown term in moment" },
    { Error::QUOTED, "Must be double-quoted" },
    { Error::LIST, "Unknown value in list" }
  } );
  
  //! error handler
  template< Error key >
  static void handleError(const std::string& value) {
    const auto& msg = err_msg.find(key);
    if (msg != err_msg.end()) {
      if (!value.empty()) {
        Throw(ExceptType::FATAL,
              "Error while parsing '" + value + "'. " + msg->second + ".");
      } else {
        Throw(ExceptType::FATAL,
              "Error while parsing. " + msg->second + ".");
      }
    } else {
      Throw(ExceptType::FATAL, "Unknown input deck parser error.");
    }
  }

  // Common actions

  //! error dispatch
  template< class Stack, Error key >
  struct error : action_base< error<Stack,key> > {
    static void apply(const std::string& value, Stack& stack) {
      handleError<key>(value);
      IGNORE(stack);    // suppress compiler warning: parameter never referenced
    }
  };

  // Common grammar

  //! read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  // match unknown keyword and handle error
  template< class Stack, Error key >
  struct unknown :
         pad< ifapply< trim<any, space>, error<Stack,key> >, blank, space > {};

  //! read 'token' padded by blank at left and space at right
  template< class token >
  struct read :
         pad< trim<token, space>, blank, space > {};

  //! parse input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions... >, blank, space > {};

  //! comment: start with '#' until eol
  struct comment :
         pad< trim<one<'#'>,eol>, blank, eol> {};

  //! number: optional sign followed by digits
  struct number :
         seq< opt< sor<one<'+'>, one<'-'>> >, digit> {};

  //! plow through 'tokens' until 'endkeyword'
  template< class Stack, typename endkeyword, typename... tokens >
  struct block :
         until< read<endkeyword>,
                sor<comment, tokens..., unknown<Stack,Error::KEYWORD>> > {};

  //! plow through list of values between keywords 'key' and "end", calling
  //! 'insert' for each if matches and allowing comments between values
  template< class Stack, typename endkeyword, class key, class insert,
            class value = number >
  struct list :
         ifmust< read<key>,
                 until< read<endkeyword>,
                        sor<comment,
                            parse<value,insert>,
                            unknown<Stack,Error::LIST>> > > {};

  //! scan string between characters 'lbound' and 'rbound' and if matches apply
  //! action 'insert'
  template< class Stack, class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         ifmust< one<lbound>,
                 ifapply< sor<trim<not_one<lbound>, one<rbound>>,
                              unknown<Stack,Error::QUOTED>>, insert >,
                 one<rbound>
               > {};

  //! process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<keywords,insert> > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class Stack, class keyword, class insert, char lbound='"',
            char rbound='"' >
  struct process_quoted :
         ifmust< read<keyword>,
                 sor< quoted<Stack,insert,lbound,rbound>,
                      unknown<Stack,Error::QUOTED>> > {};

} // grm::
} // quinoa::

#endif // QuinoaInputDeckGrammar_h
