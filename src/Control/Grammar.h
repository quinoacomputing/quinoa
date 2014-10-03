//******************************************************************************
/*!
  \file      src/Control/Grammar.h
  \author    J. Bakosi
  \date      Thu 02 Oct 2014 08:42:15 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Common of grammars
  \details   Common of grammars
*/
//******************************************************************************
#ifndef Grammar_h
#define Grammar_h

#include <sstream>

#include <Exception.h>
#include <tkTags.h>

namespace tk {
//! Grammar definition: state, actions, grammar
namespace grm {

  //! Parser's printer: this should be defined once per library in global-scope
  //! by a parser. It is defined in Control/[executable]/CmdLine/Parser.C, since
  //! every executable has at least a command line parser.
  extern Print g_print;

  // Common auxiliary functions

  //! C-style enum indicating warning or error (used as template argument)
  enum MsgType { ERROR=0, WARNING };

  //! Parser error types
  enum class MsgKey : uint8_t { KEYWORD,
                                MOMENT,
                                QUOTED,
                                LIST,
                                ALIAS,
                                MISSING,
                                UNSUPPORTED,
                                NOOPTION,
                                NOTSELECTED,
                                EXISTS,
                                NODEPVAR,
                                NOTALPHA,
                                NOTERMS,
                                NOSAMPLES,
                                INVALIDSAMPLESPACE,
                                MALFORMEDSAMPLE,
                                INVALIDBINSIZE,
                                INVALIDEXTENT,
                                EXTENTLOWER,
                                NOBINS,
                                ZEROBINSIZE,
                                MAXSAMPLES,
                                MAXBINSIZES,
                                MAXEXTENTS,
                                BINSIZES,
                                PDF,
                                PDFEXISTS,
                                BADPRECISION,
                                PRECISIONBOUNDS,
                                CHARMARG };

  //! Associate parser errors to error messages
  static const std::map< MsgKey, std::string > message( {
    { MsgKey::KEYWORD, "Unknown keyword." },
    { MsgKey::MOMENT, "Unknown term in moment." },
    { MsgKey::QUOTED, "Must be double-quoted." },
    { MsgKey::LIST, "Unknown value in list." },
    { MsgKey::ALIAS, "Alias keyword too long. Use either a full-length keyword "
      "with double-hyphens, e.g., --keyword, or its alias, a single character, "
      "with a single hyphen, e.g., -k." },
    { MsgKey::MISSING, "Required field missing." },
    { MsgKey::UNSUPPORTED, "Option not supported." },
    { MsgKey::NOOPTION, "Option does not exist." },
    { MsgKey::NOTSELECTED, "Option is not among the selected ones. The keyword "
      "here is appropriate, but in order to use this keyword in this context, "
      "the option must be selected upstream." },
    { MsgKey::EXISTS, "Dependent variable already used." },
    { MsgKey::NODEPVAR, "Dependent variable not selected. To request a "
      "statistic or PDF involving this variable, an equation must be specified "
      "upstream in the input deck assigning this dependent variable to an "
      "equation to be integrated using the depvar keyword." },
    { MsgKey::NOTALPHA, "Variable not alphanumeric." },
    { MsgKey::NOTERMS, "Statistic requires at least one variable." },
    { MsgKey::NOSAMPLES, "PDF requires at least one sample space variable." },
    { MsgKey::INVALIDSAMPLESPACE, "PDF sample space specification incorrect. A "
      "non-empty list of sample space variables, must be followed by a "
      "colon, followed by a non-empty list of bin sizes (reals numbers), e.g., "
      "\"(x y : 0.1 0.2)\"" },
    { MsgKey::MALFORMEDSAMPLE, "A PDF sample space variable must be a single "
      "upper or lowercase letter optionally followed by a single digit. "
      "Multiple variables, specifying a multi-dimensional sample space, must "
      "be separated by white spaces." },
    { MsgKey::INVALIDBINSIZE, "PDF sample space bin size(s) specification "
      "incorrect. A non-empty list of sample space variables, must be followed "
      "by a colon, followed by a non-empty list of bin sizes (real numbers), "
      "e.g., \"(x y : 0.1 0.2)\"" },
    { MsgKey::INVALIDEXTENT, "PDF sample space extents specification "
      "incorrect. The semi-colon following the list of bin sizes, must be "
      "followed by a non-empty list of extents (real numbers), e.g., \"(x y : "
      "0.1 0.2 ; 0.0 1.0 0.2 0.9)\". The number of real numbers representing "
      "the sample space extents must be exactly twice the number of sample "
      "space dimensions, i.e., in this 2D example 4 (2 pairs)." },
    { MsgKey::EXTENTLOWER, "PDF sample space extents must be a pair of a "
      "smaller and a larger numerical value, in that order." },
    { MsgKey::NOBINS, "Need at least one sample space bin size, followed by a "
      "colon, in a PDF specification." },
    { MsgKey::ZEROBINSIZE, "Sample space bin size must be a real number and "
      "greater than zero." },
    { MsgKey::MAXSAMPLES, "The maximum number of sample space variables for a "
      "joint PDF is 3." },
    { MsgKey::MAXBINSIZES, "The maximum number of bins sizes for a joint PDF "
      "is 3."},
    { MsgKey::MAXEXTENTS, "The maximum number of optional sample space extents "
      "for a joint PDF is 3 pairs."},
    { MsgKey::BINSIZES, "The number of sample space variables for a PDF must "
      "equal the number of bin sizes given." },
    { MsgKey::PDF, "Syntax error while parsing PDF specification." },
    { MsgKey::PDFEXISTS, "PDF already exists. PDF identifiers must be unique."},
    { MsgKey::BADPRECISION, "Precision specification invalid. It should be a "
      "positive integer or the word \'max\', selecting the maximum number of "
      "digits for the underyling floating point type."},
    { MsgKey::PRECISIONBOUNDS, "Precision specification out of bounds. It "
      "should be a positive integer between 1 and the maximum number of digits "
      "for the underyling floating point type on the machine. (Set \'max\' for "
      "the maximum.)"},
    { MsgKey::CHARMARG, "Arguments starting with '+' are assumed to be inteded "
      "for the Charm++ runtime system. Did you forget to prefix the command "
      "line with charmrun? If this warning persists even after running with "
      "charmrun, then Charm++ does not understand it either. See the Charm++ "
      "manual at http://charm.cs.illinois.edu/manuals/html/charm++/"
      "manual.html." }
  } );

  //! parser error and warning message handler
  template< class Stack, MsgType type, MsgKey key >
  static void Message( Stack& stack, const std::string& value ) {
    const auto& msg = message.find(key);
    if (msg != message.end()) {
      std::stringstream ss;
      const std::string typestr( type == MsgType::ERROR ? "Error" : "Warning" );
      if (!value.empty()) {
        ss << typestr << " while parsing '" << value << "' at "
           << stack.location() << ". " << msg->second;
      } else {
        ss << typestr << " while parsing at " << stack.location() << ". "
           << msg->second;
      }
      stack.template push_back< tag::error >( ss.str() );
    } else {
      stack.template push_back< tag::error >
        ( std::string("Unknown parser ") + 
          (type == MsgType::ERROR ? "error" : "warning" ) +
          " with no location information." );
    }
  }

  //! put option in state at position given by tags
  template< class Stack, class Option, class DefaultStack, class... tags >
  static void store_option( Stack& stack,
                            const std::string& value,
                            const DefaultStack& defaults ) {
    Option opt;
    //! Emit warning on overwrite
    if (stack.template get< tags... >() != defaults.template get< tags... >()) {
      g_print << "\n>>> WARNING: Multiple definitions for '"
              << opt.group() << "' option. Overwriting '"
              << opt.name( stack.template get< tags... >() ) << "' with '"
              << opt.name( opt.value( value ) ) << "'.\n\n";
    }
    if (opt.exist(value))
      stack.template set< tags... >( opt.value( value ) );
    else
      Message< Stack, ERROR, MsgKey::NOOPTION >( stack, value );
  }

  // Common PEGTL actions

  //! error message dispatch
  template< class Stack, MsgKey key >
  struct error : pegtl::action_base< error<Stack,key> > {
    static void apply(const std::string& value, Stack& stack) {
      Message< Stack, ERROR, key >( stack, value );
    }
  };

  //! warning message dispatch
  template< class Stack, MsgKey key >
  struct warning : pegtl::action_base< warning<Stack,key> > {
    static void apply(const std::string& value, Stack& stack) {
      Message< Stack, WARNING, key >( stack, value );
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

  //! convert and push back value to vector of back of vector in state at
  //! position given by tags
  template< class Stack, typename tag, typename...tags >
  struct Store_back_back :
  pegtl::action_base< Store_back_back< Stack, tag, tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      stack.template store_back_back< tag, tags... >( value );
    }
  };

  //! put true in switch in state at position given by tags
  template< class Stack, typename tag, typename... tags >
  struct Store_switch : pegtl::action_base< Store_switch<Stack,tag,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.template set<tag,tags...>(true);
    }
  };

  //! push back option in state at position given by tags
  template< class Stack, class Option, typename tag, typename... tags >
  struct store_back_option :
  pegtl::action_base< store_back_option<Stack, Option, tag, tags...> > {
    static void apply( const std::string& value, Stack& stack ) {
      Option opt;
      if (opt.exist(value))
        stack.template push_back<tag,tags...>( opt.value( value ) );
      else
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, value );
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
      stack.template insert_field<key_type, field, tag, tags...>( key, value );
    }
  };

  //! convert and insert option value to map at position given by tags
  template< class Stack, class Option, typename field, typename sel,
            typename vec, typename tag, typename... tags >
  struct insert_option :
  pegtl::action_base< insert_option< Stack, Option, field, sel, vec, tag,
                                     tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      stack.template
        insert_opt< key_type, field, typename Option::EnumType, tag, tags... >
                  ( key, Option().value(value) );
    }
  };

  // Common grammar

  //! read 'token' until 'erased' trimming, i.e., not consuming, 'erased'
  template< class token, class erased >
  struct trim :
         pegtl::seq< token, pegtl::until< pegtl::at< erased > > > {};

  // match unknown keyword and handle error
  template< class Stack, MsgKey key, template< class, MsgKey > class msg >
  struct unknown :
         pegtl::pad< pegtl::ifapply< trim< pegtl::any, pegtl::space >,
                                     msg< Stack, key > >,
                     pegtl::blank,
                     pegtl::space > {};

  //! match alias cmdline keyword
  template< class Stack, class keyword >
  struct alias :
         pegtl::seq<
           pegtl::one<'-'>,
           typename keyword::pegtl_alias,
           pegtl::sor< pegtl::space,
                       pegtl::apply< error< Stack, MsgKey::ALIAS > > > > {};

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
  //! 'keywords', apply 'actions'; as opposed to scan_until this rule, allows
  //! multiple actions
  template< class keywords, class... actions >
  struct scan :
         pegtl::pad< pegtl::ifapply< trim< keywords, pegtl::space >,
                                     actions... >,
                     pegtl::blank,
                     pegtl::space > {};

  //! scan input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'action'; using additional custom end rule, as opposed
  //! to scan, this rule allows an additional end-rule until which parsing is
  //! continued, the additional custom end-rule is OR'd to pegtl::space
  template< class keywords, class action, class end = pegtl::space >
  struct scan_until :
         pegtl::pad< pegtl::ifapply< trim< keywords,
                                           pegtl::sor< pegtl::space, end > >,
                                     action >,
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

  //! plow through 'tokens' until 'endkeyword'
  template< class Stack, class endkeyword, typename... tokens >
  struct block :
         pegtl::until<
           readkw< typename endkeyword::pegtl_string >,
           pegtl::sor< comment,
                       tokens...,
                       unknown< Stack, MsgKey::KEYWORD, error > > > {};

  //! plow through vector of values between keywords 'key' and 'endkeyword',
  //! calling 'insert' for each if matches and allowing comments between values
  template< class Stack, class key, class insert, class endkeyword,
            class starter, class value = number >
  struct vector :
         pegtl::ifmust< readkw< key >,
                        starter,
                        pegtl::until<
                          readkw< typename endkeyword::pegtl_string >,
                          pegtl::sor<
                            comment,
                            scan< value, insert >,
                            unknown< Stack, MsgKey::LIST, error > > > > {};

  //! scan string between characters 'lbound' and 'rbound' and if matches apply
  //! action 'insert'
  template< class Stack, class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         pegtl::ifmust< pegtl::one< lbound >,
                        pegtl::ifapply<
                          pegtl::sor< trim< pegtl::not_one< lbound >,
                                            pegtl::one< rbound > >,
                                      unknown< Stack, MsgKey::QUOTED, error > >,
                        insert >,
                        pegtl::one< rbound > > {};

  //! process 'keyword' and call its 'insert' action if matches 'kw_type'
  template< class Stack, class keyword, class insert,
            class kw_type = pegtl::digit >
  struct process :
         pegtl::ifmust< readkw< keyword >,
                        scan< pegtl::sor<
                                kw_type,
                                pegtl::apply< error< Stack,
                                                     MsgKey::MISSING > > >,
                              insert > > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class Stack, class keyword, class insert, char lbound='"',
            char rbound='"' >
  struct process_quoted :
         pegtl::ifmust< readkw< keyword >,
                        pegtl::sor<
                          quoted< Stack, insert, lbound, rbound >,
                          unknown< Stack, MsgKey::QUOTED, error > > > {};

  //! process command line 'keyword' and call its 'insert' action if matches
  //! 'kw_type'
  template< class Stack, class keyword, class insert,
            class kw_type = pegtl::any >
  struct process_cmd :
         pegtl::ifmust<
           readcmd< Stack, keyword >,
           scan< pegtl::sor< kw_type,
                             pegtl::apply< error< Stack, MsgKey::MISSING > > >,
                 insert > > {};

  //! process command line switch 'keyword'
  template< class Stack, class keyword, typename tag, typename... tags >
  struct process_cmd_switch :
         pegtl::ifmust<
           readcmd< Stack, keyword >,
           pegtl::apply< Store_switch< Stack, tag, tags... > > > {};

  //! read_file entry point: parse 'keywords' and 'ignore' until eof
  template< class Stack, typename keywords, typename ignore >
  struct read_file :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         ignore,
                         unknown< Stack, MsgKey::KEYWORD, error > > > {};

  //! process but ignore Charm++'s charmrun arguments starting with '+'
  template< class Stack >
  struct charmarg :
         pegtl::seq< pegtl::one<'+'>,
                     unknown< Stack, MsgKey::CHARMARG, warning > > {};

  //! read_string entry point: parse 'keywords' until end of string
  template< class Stack, typename keywords >
  struct read_string :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         charmarg< Stack >,
                         unknown< Stack, MsgKey::KEYWORD, error > > > {};

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
