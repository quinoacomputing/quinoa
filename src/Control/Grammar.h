//******************************************************************************
/*!
  \file      src/Control/Grammar.h
  \author    J. Bakosi
  \date      Wed Feb 19 15:29:20 2014
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

  // Common auxiliary functions

  //! Parser error types
  enum class Error : uint8_t { KEYWORD,
                               MOMENT,
                               QUOTED,
                               LIST,
                               ALIAS,
                               MISSING,
                               UNSUPPORTED,
                               NOOPTION,
                               NOTSELECTED,
                               EXISTS,
                               NOTALPHA };

  static const std::map< Error, std::string > err_msg( {
    { Error::KEYWORD, "Unknown keyword" },
    { Error::MOMENT, "Unknown term in moment" },
    { Error::QUOTED, "Must be double-quoted" },
    { Error::LIST, "Unknown value in list" },
    { Error::ALIAS, "Alias keyword too long" },
    { Error::MISSING, "Required field missing" },
    { Error::UNSUPPORTED, "Option not supported" },
    { Error::NOOPTION, "Option does not exist" },
    { Error::NOTSELECTED, "Option is not among the selected ones" },
    { Error::EXISTS, "Dependent variable already used" },
    { Error::NOTALPHA, "Variable not alphanumeric" }
  } );
  
  //! parser error handler
  template< class Stack, Error key >
  static void handleError( const Stack& stack, const std::string& value ) {
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

  //! put option in state at position given by tags
  template< class Stack, class OptionType, class DefaultStack, class... tags >
  static void store_option( Stack& stack,
                            const std::string& value,
                            const DefaultStack& defaults ) {
    tk::Option< OptionType > opt;
    //! Emit warning on overwrite
    if (stack.template get< tags... >() != defaults.template get< tags... >()) {
      std::cout << "\n>>> PARSER WARNING: Multiple definitions for '"
                << opt.group() << "' option. Overwriting '"
                << opt.name( stack.template get< tags... >() ) << "' with '"
                << opt.name( opt.value( value ) ) << "'.\n\n";
    }
    stack.template set< tags... >( opt.value( value ) );
  }

  // Common PEGTL actions

  //! error dispatch
  template< class Stack, Error key >
  struct error : pegtl::action_base< error<Stack,key> > {
    static void apply(const std::string& value, Stack& stack) {
      handleError< Stack, key >( stack, value );
    }
  };

  //! put value in state at position given by tags without conversion
  template<class Stack, typename tag, typename... tags >
  struct Set : pegtl::action_base< Set<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template set<tag,tags...>(value);
    }
  };

  //! put value in state at position given by tags with conversion, see Control
  template< class Stack, typename tag, typename... tags >
  struct Store : pegtl::action_base< Store<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template store<tag,tags...>(value);
    }
  };

  //! convert and push back value to vector in state at position given by tags
  template< class Stack, typename tag, typename...tags >
  struct Store_back : pegtl::action_base< Store_back<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template store_back<tag,tags...>(value);
    }
  };

  //! push back option in state at position given by tags
  template< class Stack, class OptionType, typename tag, typename... tags >
  struct store_back_option :
  pegtl::action_base< store_back_option<Stack, OptionType, tag, tags...> > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::Option< OptionType > opt;
      if (opt.exist( value ))
        stack.template push_back<tag,tags...>( opt.value( value ) );
      else
        handleError< Stack, Error::NOOPTION >( stack, value );
    }
  };

  //! convert and insert value to map at position given by tags
  template< class Stack, typename field, typename sel, typename vec,
            typename tag, typename...tags >
  struct Insert_field :
  pegtl::action_base< Insert_field< Stack, field, sel, vec, tag, tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      stack.template insert_field< key_type, field, tag, tags... >( key, value );
    }
  };

  //! convert and insert option value to map at position given by tags
  template< class Stack, class OptionType, typename field, typename sel,
            typename vec, typename tag, typename... tags >
  struct insert_option :
  pegtl::action_base< insert_option< Stack, OptionType, field, sel, vec, tag,
                                     tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::Option< OptionType > opt;
      using EnumType = typename OptionType::EnumType;
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      stack.template insert_opt< key_type, field, EnumType, tag, tags... >
                               ( key, opt.value(value) );
    }
  };

  // Common grammar

  //! read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         pegtl::seq< token, pegtl::until< pegtl::at< erased > > > {};

  // match unknown keyword and handle error
  template< class Stack, Error key >
  struct unknown :
         pegtl::pad< pegtl::ifapply< trim< pegtl::any, pegtl::space >,
                                     error< Stack, key > >,
                     pegtl::blank,
                     pegtl::space > {};

  //! match alias cmdline keyword
  template< class Stack, class keyword >
  struct alias :
         pegtl::seq< pegtl::one<'-'>,
                     typename keyword::pegtl_alias,
                     pegtl::sor< pegtl::space,
                                 unknown< Stack, Error::ALIAS > > > {};

  //! match verbose cmdline keyword
  template< class keyword >
  struct verbose :
         pegtl::seq< pegtl::string<'-','-'>,
                     typename keyword::pegtl_string,
                     pegtl::space > {};

  //! read keyword 'token' padded by blank at left and space at right
  template< class token >
  struct readkw :
         pegtl::pad< trim< token, pegtl::space >,
                     pegtl::blank,
                     pegtl::space > {};

  //! read command line 'keyword' in either verbose or alias form
  template< class Stack, class keyword >
  struct readcmd :
         pegtl::sor< verbose< keyword >, alias< Stack, keyword > > {};

  //! scan input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct scan :
         pegtl::pad< pegtl::ifapply< trim< keywords, pegtl::space >,
                                     actions... >,
                     pegtl::blank,
                     pegtl::space > {};

  //! comment: start with '#' until eol
  struct comment :
         pegtl::pad< trim< pegtl::one<'#'>, pegtl::eol >,
                     pegtl::blank,
                     pegtl::eol > {};

  //! number: optional sign followed by digits
  struct number :
         pegtl::seq< pegtl::opt< pegtl::sor< pegtl::one<'+'>,
                                             pegtl::one<'-'> > >,
                     pegtl::digit > {};

  //! plow through 'tokens' until kw::end
  template< class Stack, typename... tokens >
  struct block :
         pegtl::until< readkw< kw::end::pegtl_string >,
                       pegtl::sor< comment,
                                   tokens...,
                                   unknown< Stack, Error::KEYWORD > > > {};

  //! plow through vector of values between keywords 'key' and kw::end, calling
  //! 'insert' for each if matches and allowing comments between values
  template< class Stack, class key, class insert, class value = number >
  struct vector :
         pegtl::ifmust< readkw< key >,
                        pegtl::until< readkw< kw::end::pegtl_string >,
                                      pegtl::sor< comment,
                                                  scan< value, insert >,
                                                  unknown<Stack,
                                                          Error::LIST> > > > {};

  //! scan string between characters 'lbound' and 'rbound' and if matches apply
  //! action 'insert'
  template< class Stack, class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         pegtl::ifmust< pegtl::one< lbound >,
                        pegtl::ifapply<
                          pegtl::sor< trim< pegtl::not_one< lbound >,
                                            pegtl::one< rbound > >,
                                      unknown< Stack, Error::QUOTED > >,
                        insert >,
                        pegtl::one< rbound > > {};

  //! process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class Stack, class keyword, class insert,
            class keywords = pegtl::digit >
  struct process :
         pegtl::ifmust< readkw< keyword >,
                        scan< pegtl::sor<
                                keywords,
                                pegtl::apply< error< Stack,
                                                     Error::MISSING > > >,
                              insert > > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class Stack, class keyword, class insert, char lbound='"',
            char rbound='"' >
  struct process_quoted :
         pegtl::ifmust< readkw< keyword >,
                        pegtl::sor< quoted< Stack, insert, lbound, rbound >,
                                    unknown< Stack, Error::QUOTED > > > {};

  //! process command line 'keyword' and call its 'insert' action if matches
  //! 'keywords'
  template< class Stack, class keyword, class insert,
            class keywords = pegtl::any >
  struct process_cmd :
         pegtl::ifmust< readcmd< Stack, keyword >,
                        scan< pegtl::sor<
                                keywords,
                                pegtl::apply< error< Stack,
                                                     Error::MISSING > > >,
                              insert > > {};

  //! read_file entry point: parse 'keywords' and 'ignore' until eof
  template< class Stack, typename keywords, typename ignore >
  struct read_file :
         pegtl::until< pegtl::eof,
                       pegtl::sor< keywords,
                                   ignore,
                                   unknown< Stack, Error::KEYWORD > > > {};

  //! read_string entry point: parse 'keywords' until end of string
  template< class Stack, typename keywords >
  struct read_string :
         pegtl::until< pegtl::eof,
                       pegtl::sor< keywords,
                                   unknown< Stack, Error::KEYWORD > > > {};

  //! insert RNG parameter
  template< typename Stack, typename keyword, typename option, typename field,
            typename sel, typename vec, typename... tags >
  struct rng_option :
         process< Stack,
                  typename keyword::pegtl_string,
                  insert_option< Stack,
                                 option,
                                 field,
                                 sel, vec, tags... >,
                  pegtl::alpha > {};

} // grm::
} // tk::

#endif // Grammar_h
