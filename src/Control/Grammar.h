//******************************************************************************
/*!
  \file      src/Control/Grammar.h
  \author    J. Bakosi
  \date      Tue Oct 29 07:27:24 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Common of grammars
  \details   Common of grammars
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <sstream>

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
    { Error::MISSING, "Required field missing" }
  } );
  
  //! error handler
  template< class Stack, Error key >
  static void handleError(const Stack& stack, const std::string& value) {
    const auto& msg = err_msg.find(key);
    if (msg != err_msg.end()) {
      if (!value.empty()) {
        std::stringstream ss;
        ss << "Error while parsing '" << value << "' at " << stack.location()
           << ". " << msg->second << ".";
        Throw(ExceptType::FATAL, ss.str());
      } else {
        std::stringstream ss;
        ss << "Error while parsing at " << stack.location() << ". "
           << msg->second << ".";
        Throw(ExceptType::FATAL, ss.str());
      }
    } else {
      Throw(ExceptType::FATAL, "Unknown parser error.");
    }
  }

  // Common actions

  //! error dispatch
  template< class Stack, Error key >
  struct error : action_base< error<Stack,key> > {
    static void apply(const std::string& value, Stack& stack) {
      handleError<Stack,key>(stack,value);
    }
  };

  //! put value in state at position given by tags without conversion
  template<class Stack, typename tag, typename... tags >
  struct Set : action_base< Set<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template set<tag,tags...>(value);
    }
  };

  //! put value in state at position given by tags with conversion, see Control
  template< class Stack, typename tag, typename... tags >
  struct Store : action_base< Store<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template store<tag,tags...>(value);
    }
  };

  //! convert and push back value to vector in state at position given by tags
  template< class Stack, typename tag, typename...tags >
  struct Store_back : action_base< Store_back<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template store_back<tag,tags...>(value);
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

  //! match alias cmdline keyword
  template< class Stack, class keyword >
  struct alias :
         seq< one<'-'>,
              typename keyword::pegtl_alias,
              sor<space, unknown<Stack,Error::ALIAS>> > {};

  //! match verbose cmdline keyword
  template< class keyword >
  struct verbose :
         seq< string<'-','-'>, typename keyword::pegtl_string, space > {};

  //! read keyword 'token' padded by blank at left and space at right
  template< class token >
  struct readkw :
         pad< trim<token, space>, blank, space > {};

  //! read command line 'keyword' in either verbose or alias form
  template< class Stack, class keyword >
  struct readcmd :
         sor< verbose<keyword>, alias<Stack,keyword> > {};

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
  template< class Stack, class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< readkw<keyword>,
                 scan< sor<keywords, apply<error<Stack,Error::MISSING>>>,
                       insert> > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class Stack, class keyword, class insert, char lbound='"',
            char rbound='"' >
  struct process_quoted :
         ifmust< readkw<keyword>,
                 sor< quoted<Stack,insert,lbound,rbound>,
                      unknown<Stack,Error::QUOTED>> > {};

  //! process command line 'keyword' and call its 'insert' action if matches
  //! 'keywords'
  template< class Stack, class keyword, class insert, class keywords = any >
  struct process_cmd :
         ifmust< readcmd<Stack,keyword>,
                 scan< sor<keywords, apply<error<Stack,Error::MISSING>>>,
                       insert> > {};

} // grm::
} // tk::

#endif // QuinoaInputDeckGrammar_h
