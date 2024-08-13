// *****************************************************************************
/*!
  \file      src/Control/CommonGrammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Generic, low-level grammar, re-used by specific grammars
  \details   Generic, low-level grammar. We use the Parsing Expression Grammar
    Template Library (PEGTL) to create the grammar and the associated parser.
*/
// *****************************************************************************
#ifndef CommonGrammar_h
#define CommonGrammar_h

#include <type_traits>
#include <sstream>

#include <brigand/algorithms/for_each.hpp>
#include <brigand/functions/logical/or.hpp>
#include <brigand/sequences/has_key.hpp>

#include "If.hpp"
#include "Exception.hpp"
#include "Tags.hpp"
#include "StatCtr.hpp"
#include "Options/TxtFloatFormat.hpp"
#include "Options/Error.hpp"

namespace tk {
//! Toolkit general purpose grammar definition
namespace grm {

  using namespace tao;

  using ncomp_t = kw::ncomp::info::expect::type;

  //! Parser's printer: this should be defined once per library in global-scope
  //! (still in namespace, of course) by a parser. It is defined in
  //! Control/[executable]/CmdLine/Parser.C, since every executable has at least
  //! a command line parser.
  extern Print g_print;

  // Common auxiliary functions (reused by multiple grammars)

  //! C-style enum indicating warning or error (used as template argument)
  enum MsgType { ERROR=0, WARNING };

  //! Parser error types
  enum class MsgKey : uint8_t {
    KEYWORD,            //!< Unknown keyword
    MOMENT,             //!< Unknown Term in a Moment
    QUOTED,             //!< String must be double-quoted
    LIST,               //!< Unknown value in list
    ALIAS,              //!< Alias keyword too long
    MISSING,            //!< Required field missing
    PREMATURE,          //!< Premature end of line
    UNSUPPORTED,        //!< Option not supported
    NOOPTION,           //!< Option does not exist
    NOINIT,             //!< No (or too many) initialization policy selected
    NOPROBLEM,          //!< No test problem type selected
    NOCOEFF,            //!< No coefficients policy selected
    NOTSELECTED,        //!< Option not selected upstream
    NOSOLVE,            //!< Dependent variable to solve for has not been spec'd
    NOTALPHA,           //!< Variable must be alphanumeric
    POINTEXISTS,        //!< Point identifier already defined
    BADPRECISION,       //!< Floating point precision specification incorrect
    BOUNDS,             //!< Specified value out of bounds
    PRECISIONBOUNDS,    //!< Floating point precision spec out of bounds
    UNFINISHED,         //!< Unfinished block
    CHARMARG,           //!< Argument inteded for the Charm++ runtime system
    OPTIONAL };         //!< Message key used to indicate of something optional

  //! Associate parser errors to error messages
  static const std::map< MsgKey, std::string > message{
    { MsgKey::KEYWORD, "Unknown keyword or keyword unrecognized in this "
      "block." },
    { MsgKey::MOMENT, "Unknown term in moment." },
    { MsgKey::QUOTED, "Must be double-quoted." },
    { MsgKey::LIST, "Unknown value in list." },
    { MsgKey::ALIAS, "Alias keyword too long. Use either a full-length keyword "
      "with double-hyphens, e.g., --keyword, or its alias, a single character, "
      "with a single hyphen, e.g., -k." },
    { MsgKey::MISSING, "Required field missing." },
    { MsgKey::PREMATURE, "Premature end of line." },
    { MsgKey::UNSUPPORTED, "Option not supported." },
    { MsgKey::NOOPTION, "Option does not exist." },
    { MsgKey::NOTSELECTED, "Option is not among the selected ones. The keyword "
      "here is appropriate, but in order to use this keyword in this context, "
      "the option must be selected upstream." },
    { MsgKey::NOTALPHA, "Variable not alphanumeric." },
    { MsgKey::NOSOLVE, "Dependent variable to solve for not specified within "
      "the block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'solve' to specify the type of the dependent "
      "variable to solve for." },
    { MsgKey::NOINIT, "No (or too many) initialization policy (or policies) "
      "has been specified within the block preceding this position. An "
      "initialization policy (and only one) is mandatory for the preceding "
      "block. Use the keyword 'init' to specify an initialization policy." },
    { MsgKey::NOPROBLEM, "No test problem has been specified within the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'problem' to specify a test problem." },
    { MsgKey::NOCOEFF, "No coefficients policy has been specified within the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'coeff' to specify an coefficients policy." },
    { MsgKey::POINTEXISTS, "Point already exists. Point identifiers must be "
      "unique."},
    { MsgKey::BADPRECISION, "Precision specification invalid. It should be a "
      "positive integer or the word \'max\', selecting the maximum number of "
      "digits for the underyling floating point type."},
    { MsgKey::BOUNDS, "Specified value out of bounds. For the bounds on a "
      "keyword, run '<executable> -H <keyword>'."},
    { MsgKey::PRECISIONBOUNDS, "Precision specification out of bounds. It "
      "should be a positive integer between 1 and the maximum number of digits "
      "for the underyling floating point type on the machine. (Set \'max\' for "
      "the maximum.)"},
    { MsgKey::UNFINISHED, "Block started but not finished by the 'end' "
      "keyword." },
    { MsgKey::CHARMARG, "Arguments starting with '+' are assumed to be inteded "
      "for the Charm++ runtime system. Did you forget to prefix the command "
      "line with charmrun? If this warning persists even after running with "
      "charmrun, then Charm++ does not understand it either. See the Charm++ "
      "manual at http://charm.cs.illinois.edu/manuals/html/charm++/"
      "manual.html." },
    { MsgKey::OPTIONAL, "This is not really an error message and thus it "
      "should not be used as one. But its key can be used to indicate "
      "something optional (which is not an error), which in some situations is "
      "not optional (which is an error)." }
  };

  //! Parser error and warning message handler.
  //! \details This function is used to associated and dispatch an error or a
  //!   warning during parsing. After finding the error message corresponding to
  //!   a key, it pushes back the message to a std::vector of std::string, which
  //!   then will be diagnosed later by tk::FileParser::diagnostics. The
  //!   template arguments define (1) the grammar stack (Stack, a tagged tuple)
  //!   to operate on, (2) the message type (error or warning), and (3) the
  //!   message key used to look up the error message associated with the key.
  //! \param[in,out] stack Grammar stack (a tagged tuple) to operate on
  //! \param[in] in Last parsed PEGTL input token (can be empty, depending on
  //!   what context this function gets called.
  template< class Stack, MsgType type, MsgKey key, class Input >
  static void Message( Stack& stack, const Input& in ) {
    const auto& msg = message.find(key);
    if (msg != message.end()) {
      std::stringstream ss;
      const std::string typestr( type == MsgType::ERROR ? "Error" : "Warning" );
      auto pos = in.position();
      if (!in.empty()) {
        ss << typestr << " while parsing '" << in.string() << "' at "
           << pos.line << ',' << pos.byte_in_line << ". " << msg->second;
      } else {
        ss << typestr << " while parsing at " << pos.line << ','
           << pos.byte_in_line << ". " << msg->second;
      }
      stack.template get< tag::error >().push_back( ss.str() );
    } else {
      stack.template get< tag::error >().push_back(
        std::string("Unknown parser ") +
        (type == MsgType::ERROR ? "error" : "warning" ) +
        " with no location information." );
    }
  }

  // Common PEGTL actions (PEGTL actions reused by multiple grammars)

  //! PEGTL action base: do nothing by default
  //! \details This base is specialized to different actions duing parsing.
  template< typename Rule >
  struct action : pegtl::nothing< Rule > {};

  //! Helper for calling action::apply for multiple actions
  template< typename... As >
  struct call {
    template< typename Input, typename State >
    static void apply( const Input& in, State& state ) {
        using swallow = bool[];
        (void)swallow{ ( action< As >::apply( in, state ), true )..., true };
     }
  };

  //! Rule used to trigger action(s) for a rule
  template< class rule, class... actions >
  struct act : rule {};

  //! \details Specialization of action for act< rule, actions... >
  template< class rule, class... actions >
  struct action< act< rule, actions... > > : call< actions... > {};

  //! Rule used to trigger action
  template< MsgType, MsgKey > struct msg : pegtl::success {};
  //! Error message dispatch
  //! \details This struct and its apply function are used to dispatch a message
  //!   (e.g., error, waring) from the parser. It is simply an interface to
  //!   Message. See This struct is practically used as a functor, i.e., a
  //!   struct or class that defines the function call operator, but instead the
  //!   function call operator, PEGTL uses the apply() member function for the
  //!   actions. Thus this struct can be passed to, i.e., specialize a template,
  //!   such as tk::grm::unknown, injecting in it the desired behavior.
  template< MsgType type, MsgKey key >
  struct action< msg< type, key > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Message< Stack, type, key >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags > struct Store : pegtl::success {};
  //! Put value in state at position given by tags with conversion
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store member function of the underlying grammar
  //!    stack, tk::Control::store.
  template< typename tag, typename... tags >
  struct action< Store< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      if (!in.string().empty())
        stack.template store< tag, tags... >( in.string() );
      else
        Message< Stack, ERROR, MsgKey::MISSING >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< typename... tags >
  struct Invert_switch : pegtl::success {};
  //! Invert bool in switch at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for setting a boolean value to true in the underlying grammar
  //!    stack via the member function tk::Control::set.
  template< typename... tags >
  struct action< Invert_switch< tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template get< tags... >() = !stack.template get< tags... >();
    }
  };

  //! Rule used to trigger action
  struct helpkw : pegtl::success {};
  //! \brief Find keyword among all keywords and if found, store the keyword
  //!    and its info on which help was requested behind tag::helpkw in Stack
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper to search for a keyword in the pool of registered keywords
  //!    recognized by a grammar and store the keyword and its info on which
  //!    help was requested behind tag::helpkw. Note that this functor assumes
  //!    a specific location for the std::maps of the command-line and control
  //!    file keywords pools (behind tag::cmdinfo and tag::ctrinfo,
  //!    respectively), and for the keyword and its info on which help was
  //!    requested (behind tag::helpkw). This is the structure of CmdLine
  //!    objects, thus this functor should be called from command line parsers.
  template<>
  struct action< helpkw > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& cmdinfo = stack.template get< tag::cmdinfo >();
      const auto& ctrinfo = stack.template get< tag::ctrinfo >();
      auto it = cmdinfo.find( in.string() );
      if (it != cmdinfo.end()) {
        // store keyword and its info on which help was requested
        stack.template get< tag::helpkw >() = { it->first, it->second, true };
      } else {
        it = ctrinfo.find( in.string() );
        if (it != ctrinfo.end())
          // store keyword and its info on which help was requested
          stack.template get< tag::helpkw >() = { it->first, it->second, false };
        else
          Message< Stack, ERROR, MsgKey::KEYWORD >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  template< class keyword, typename tag, typename... tags >
  struct check_lower_bound : pegtl::success {};
  //! Check if value is larger than lower bound
  template< class keyword, typename tag, typename... tags >
  struct action< check_lower_bound< keyword, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto lower = keyword::info::expect::lower;
      auto val = stack.template get< tag, tags... >();
      if (val < lower) Message< Stack, ERROR, MsgKey::BOUNDS >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class keyword, typename tag, typename... tags >
  struct check_upper_bound : pegtl::success {};
  //! Check if value is lower than upper bound
  template< class keyword, typename tag, typename... tags >
  struct action< check_upper_bound< keyword, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto upper = keyword::info::expect::upper;
      auto val = stack.template get< tag, tags... >();
      if (val > upper) Message< Stack, ERROR, MsgKey::BOUNDS >( stack, in );
    }
  };

  // Common grammar (grammar that is reused by multiple grammars)

  //! Read 'token' until 'erased' trimming, i.e., not consuming, 'erased'
  template< class token, class erased >
  struct trim :
         pegtl::seq< token,
                     pegtl::sor<
                       pegtl::until< pegtl::at< erased > >,
                       msg< ERROR, MsgKey::PREMATURE > > > {};

  //! Match unknown keyword and handle error
  template< MsgType type, MsgKey key >
  struct unknown :
         pegtl::pad< pegtl::seq< trim< pegtl::any, pegtl::space >,
                                 msg< type, key > >,
                     pegtl::blank,
                     pegtl::space > {};

  //! Match alias cmdline keyword
  //! \details An alias command line keyword is prefixed by a single dash, '-'.
  template< class keyword >
  struct alias :
         pegtl::seq<
           pegtl::one< '-' >,
           typename keyword::info::alias::type,
           pegtl::sor< pegtl::space,
                       msg< ERROR, MsgKey::ALIAS > >  > {};

  //! Match verbose cmdline keyword
  //! \details A verbose command line keyword is prefixed by a double-dash,
  //!   '--'.
  template< class keyword >
  struct verbose :
         pegtl::seq< pegtl::string<'-','-'>,
                     typename keyword::pegtl_string,
                     pegtl::space > {};

  //! Read keyword 'token' padded by blank at left and space at right
  template< class token >
  struct readkw :
         pegtl::pad< trim< token, pegtl::space >,
                     pegtl::blank,
                     pegtl::space > {};

  //! Read command line 'keyword' in verbose form, i.e., '--keyword'
  //! \details This version is used if no alias is defined for the given keyword
  template< class keyword, typename = void >
  struct readcmd :
         verbose< keyword > {};

  //! Read command line 'keyword' in either verbose or alias form
  //! \details This version is used if an alias is defined for the given
  //!   keyword, in which case either the verbose or the alias form of the
  //!   keyword is matched, i.e., either '--keyword' or '-a', where 'a' is the
  //!   single-character alias for the longer 'keyword'. This is a partial
  //!   specialization of the simpler verbose-only readcmd, which attempts to
  //!   find the typedef 'alias' in keyword::info. If it finds it, it uses this
  //!   specialization. If it fails, it is [SFINAE]
  //!   (http://en.cppreference.com/w/cpp/language/sfinae), and falls back to
  //!   the verbose-only definition. Credit goes to David Rodriguez at
  //!   stackoverflow.com. This allows not having to change this client-code to
  //!   the keywords definitions: a keyword either defines an alias or not, and
  //!   the grammar here will do the right thing: if there is an alias, it will
  //!   build the grammar that optionally parses for it.
  //! \see tk::if_
  //! \see Control/Keyword.h and Control/Keywords.h
  //! \see http://en.cppreference.com/w/cpp/language/sfinae
  //! \see http://stackoverflow.com/a/11814074
  template< class keyword >
  struct readcmd< keyword,
                  typename if_< false, typename keyword::info::alias >::type > :
         pegtl::sor< verbose< keyword >, alias< keyword > > {};

  //! \brief Scan input padded by blank at left and space at right and if it
  //!   matches 'keyword', apply 'actions'
  //! \details As opposed to scan_until this rule, allows multiple actions
  template< class keyword, class... actions >
  struct scan :
         pegtl::pad< act< trim< keyword, pegtl::space >, actions... >,
                     pegtl::blank,
                     pegtl::space > {};

  //! Parse a number: an optional sign followed by digits
  struct number :
         pegtl::seq< pegtl::opt< pegtl::sor< pegtl::one<'+'>,
                                             pegtl::one<'-'> > >,
                     pegtl::digit > {};

  //! \brief Process command line 'keyword' and call its 'insert' action if
  //!   matches 'kw_type'
  template< template< class > class use, class keyword, class insert,
            class kw_type, class tag, class... tags >
  struct process_cmd :
         pegtl::if_must<
           readcmd< use< keyword > >,
           scan< pegtl::sor< kw_type, msg< ERROR, MsgKey::MISSING > >, insert >,
           typename std::conditional<
             tk::HasVar_expect_lower< typename keyword::info >::value,
             check_lower_bound< keyword, tag, tags... >,
             pegtl::success >::type,
           typename std::conditional<
             tk::HasVar_expect_upper< typename keyword::info >::value,
             check_upper_bound< keyword, tag, tags... >,
             pegtl::success >::type > {};

  //! Process command line switch 'keyword'
  //! \details The value of a command line switch is a boolean, i.e., it can be
  //!    either set or unset.
  template< template< class > class use, class keyword, typename tag,
            typename... tags >
  struct process_cmd_switch :
         pegtl::seq<readcmd< use<keyword> >, Invert_switch< tag, tags... >> {};

  //! Process but ignore Charm++'s charmrun arguments starting with '+'
  struct charmarg :
         pegtl::seq< pegtl::one<'+'>,
                     unknown< WARNING, MsgKey::CHARMARG > > {};

  //! Generic string parser entry point: parse 'keywords' until end of string
  template< typename keywords >
  struct read_string :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         charmarg,
                         unknown< ERROR, MsgKey::KEYWORD > > > {};

  //! \brief Ensure that a grammar only uses keywords from a pool of
  //!   pre-defined keywords
  //! \details In grammar definitions, every keyword should be wrapped around
  //!   this use template, which conditionally inherits the keyword type its
  //!   templated on (as defined in Control/Keywords.h) if the keyword is listed
  //!   in any of the pools of keywords passed in as the required pool template
  //!   argument and zero or more additional pools. If the keyword is not in any
  //!   of the pools, a compiler error is generated, since the struct cannot
  //!   inherit from base class 'char'. In that case, simply add the new keyword
  //!   into one of the pools of keywords corresponding to the given grammar.
  //!   The rationale behind this wrapper is to force the developer to maintain
  //!   the keywords pool for a grammar. The pools are brigand::set and are
  //!   used to provide help on command line arguments for a given executable.
  //!   They allow compile-time iteration with brigand::for_each or
  //!   generating a run-time std::map associating, e.g., keywords to their help
  //!   strings.
  //! \warning Note that an even more elegant solution to the problem this
  //!   wrapper is intended to solve is to use a metaprogram that collects all
  //!   occurrences of the keywords in a grammar. However, that does not seem to
  //!   be possible without extensive macro-trickery. Instead, if all keywords
  //!   in all grammar definitions are wrapped inside this use template (or one
  //!   of its specializations), we make sure that only those keywords can be
  //!   used by a grammar that are listed in the pool corresponding to a
  //!   grammar. However, this is still only a partial solution, since listing
  //!   more keywords in the pool than those used in the grammar is still
  //!   possible, which would result in those keywords included in, e.g., the
  //!   on-screen help generated even though some of the keywords may not be
  //!   implemented by the given grammar. So please don't abuse and don't list
  //!   keywords in the pool only if they are implemented in the grammar.
  //! \see For example usage see the template typedef
  //!   walker::cmd::use in Control/Walker/CmdLine/Grammar.h and its keywords
  //!   pool, walker::ctr::CmdLine::keywords, in
  //!   Control/Walker/CmdLine/CmdLine.h.
  //! \see http://en.cppreference.com/w/cpp/types/conditional
  //! TODO It still would be nice to generate a more developer-friendly
  //!    compiler error if the keyword is not in the pool.
  template< typename keyword, typename pool >
  struct use :
         std::conditional< brigand::or_<
                             brigand::has_key< pool, keyword > >::value,
                           keyword,
                           char >::type {};

} // grm::
} // tk::

#endif // CommonGrammar_h
