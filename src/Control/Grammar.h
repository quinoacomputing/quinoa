//******************************************************************************
/*!
  \file      src/Control/Grammar.h
  \author    J. Bakosi
  \date      Mon Oct  7 14:13:33 2013
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

namespace tk {
//! Grammar definition: state, actions, grammar
namespace grm {

  using namespace pegtl;

  // Common auxiliary functions

  //! Parser error types
  enum class Error : uint8_t { KEYWORD,
                               MOMENT,
                               QUOTED,
                               LIST,
                               ALIAS,
                               MISSING };

  static const std::map<Error, std::string> err_msg( {
    { Error::KEYWORD, "Unknown keyword" },
    { Error::MOMENT, "Unknown term in moment" },
    { Error::QUOTED, "Must be double-quoted" },
    { Error::LIST, "Unknown value in list" },
    { Error::ALIAS, "Alias keyword too long" },
    { Error::MISSING, "Required filename missing" }
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

  //! put value in state at position given by tags without conversion
  template<class Stack, typename tag, typename... tags >
  struct put : action_base< put<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template set<tag,tags...>(value);
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

  //! read keyword 'token' padded by blank at left and space at right
  template< class token >
  struct readkw :
         pad< trim<token, space>, blank, space > {};

  //! scan input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct scan :
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
         until< readkw<endkeyword>,
                sor<comment, tokens..., unknown<Stack,Error::KEYWORD>> > {};

  //! plow through vector of values between keywords 'key' and "end", calling
  //! 'insert' for each if matches and allowing comments between values
  template< class Stack, typename endkeyword, class key, class insert,
            class value = number >
  struct vector :
         ifmust< readkw<key>,
                 until< readkw<endkeyword>,
                        sor<comment,
                            scan<value,insert>,
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
         ifmust< readkw<keyword>, scan<keywords,insert> > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class Stack, class keyword, class insert, char lbound='"',
            char rbound='"' >
  struct process_quoted :
         ifmust< readkw<keyword>,
                 sor< quoted<Stack,insert,lbound,rbound>,
                      unknown<Stack,Error::QUOTED>> > {};

} // grm::
} // tk::

#endif // QuinoaInputDeckGrammar_h
