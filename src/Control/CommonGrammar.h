// *****************************************************************************
/*!
  \file      src/Control/CommonGrammar.h
  \author    J. Bakosi, D. Frey
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Generic, low-level grammar, re-used by specific grammars
  \details   Generic, low-level grammar. We use the Parsing Expression Grammar
    Template Library (PEGTL) to create the grammar and the associated parser.
*/
// *****************************************************************************
#ifndef CommonGrammar_h
#define CommonGrammar_h

#include <sstream>

// See documentation for tk::grm::use below for why these macros are here
#undef BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#define BOOST_MPL_CFG_NO_PREPROCESSED_HEADERS
#undef BOOST_MPL_LIMIT_METAFUNCTION_ARITY
#define BOOST_MPL_LIMIT_METAFUNCTION_ARITY 10

#include <boost/mpl/or.hpp>
#include "NoWarning/for_each.h"

#include "If.h"
#include "Exception.h"
#include "Tags.h"
#include "StatCtr.h"
#include "Options/PDFFile.h"
#include "Options/PDFPolicy.h"
#include "Options/PDFCentering.h"
#include "Options/TxtFloatFormat.h"

namespace tk {
//! Toolkit general purpose grammar definition
namespace grm {

  using ncomp_t = kw::ncomp::info::expect::type;

  //! Parser's printer: this should be defined once per library in global-scope
  //! (still in namespace, of course) by a parser. It is defined in
  //! Control/[executable]/CmdLine/Parser.C, since every executable has at least
  //! a command line parser.
  extern Print g_print;

  // Common InputDeck state

  //! Out-of-struct storage of field ID for pushing terms for statistics
  static ncomp_t field = 0;
  //! \brief Parser-lifetime storage for dependent variables selected.
  //! \details Used to track the dependent variable of differential equations
  //!   (i.e., models) assigned during parsing. It needs to be case insensitive
  //!   since we only care about whether the variable is selected or not and not
  //!   whether it denotes a full variable (upper case) or a fluctuation (lower
  //!   case). This is true for both inserting variables into the set as well as
  //!   at matching terms of products in parsing requested statistics.
  static std::set< char, tk::ctr::CaseInsensitiveCharLess > depvars;
  //! \brief Parser-lifetime storage for PDF names.
  //! \details Used to track the names registered  so that parsing new ones can
  //!    be required to be unique.
  static std::set< std::string > pdfnames;

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
    NOINIT,             //!< No initialization policy selected
    NOPROBLEM,          //!< No test problem type selected
    NOCOEFF,            //!< No coefficients policy selected
    NOTSELECTED,        //!< Option not selected upstream
    EXISTS,             //!< Variable already used
    NODEPVAR,           //!< Dependent variable has not been specified
    NOSUCHDEPVAR,       //!< Dependent variable has not been previously selected
    NOTALPHA,           //!< Variable must be alphanumeric
    NOTERMS,            //!< Statistic need a variable
    ODDSPIKES,          //!< Incomplete spikes block
    HEIGHTSPIKES,       //!< Height-sum of spikes does not add up to unity
    NODELTA,            //!< No icdelta...end block when initpolicy = jointdelta
    NOBETA,             //!< No icbeta...end block when initpolicy = jointbeta
    WRONGBETAPDF,       //!< Wrong number of parameters configuring a beta pdf
    NONCOMP,            //!< No number of components selected
    NORNG,              //!< No RNG selected
    NODT,               //!< No time-step-size policy selected
    MULDT,              //!< Multiple time-step-size policies selected
    NOSAMPLES,          //!< PDF need a variable
    INVALIDSAMPLESPACE, //!< PDF sample space specification incorrect
    MALFORMEDSAMPLE,    //!< PDF sample space variable specification incorrect
    INVALIDBINSIZE,     //!< PDF sample space bin size specification incorrect
    INVALIDEXTENT,      //!< PDF sample space extent specification incorrect
    EXTENTLOWER,        //!< PDF sample space extent-pair in non-increasing order
    NOBINS,             //!< PDF sample space bin size required
    ZEROBINSIZE,        //!< PDF sample space bin size incorrect
    MAXSAMPLES,         //!< PDF sample space dimension too large
    MAXBINSIZES,        //!< PDF sample space bin sizes too many
    MAXEXTENTS,         //!< PDF sample space extent-pairs too many
    BINSIZES,           //!< PDF sample space variables unequal to number of bins
    PDF,                //!< PDF specification syntax error
    PDFEXISTS,          //!< PDF identifier already defined
    BADPRECISION,       //!< Floating point precision specification incorrect
    PRECISIONBOUNDS,    //!< Floating point precision specification out of bounds
    UNFINISHED,         //!< Unfinished block
    VORTICAL_UNFINISHED,//!< Vortical flow problem configuration unfinished
    BC_EMPTY,           //!< Empty boundary condition block
    WRONGSIZE,          //!< Size of parameter vector incorrect
    HYDROTIMESCALES,    //!< Missing required hydrotimescales vector
    HYDROPRODUCTIONS,   //!< Missing required hydroproductions vector
    CHARMARG };         //!< Argument inteded for the Charm++ runtime system

  //! Associate parser errors to error messages
  static const std::map< MsgKey, std::string > message( {
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
    { MsgKey::EXISTS, "Dependent variable already used." },
    { MsgKey::NOSUCHDEPVAR, "Dependent variable not selected. To request a "
      "statistic or PDF involving this variable, or use this variable as a "
      "coefficients policy variable, an equation must be specified "
      "upstream in the control file assigning this variable to an "
      "equation to be integrated using the depvar keyword." },
    { MsgKey::NOTALPHA, "Variable not alphanumeric." },
    { MsgKey::HEIGHTSPIKES, "The sum of all spike heights given in the "
      "spike...end block does not add up to unity. A spike...end block "
      "must contain an even number of real numbers, where every odd one is the "
      "sample space position of a spike followed by the spike height "
      "specifying the relative probability of the spike. Since the spike "
      "heights are probabilities relative to unity, they must sum to one." },
    { MsgKey::NODEPVAR, "Dependent variable not specified within the block "
      "preceding this position. This is mandatory for the preceding block. Use "
      "the keyword 'depvar' to specify the dependent variable." },
    { MsgKey::NONCOMP, "The number of components has not been specified in the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'ncomp' to specify the number of components." },
    { MsgKey::NORNG, "The random number generator has not been specified in "
      "the block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'rng' to specify the random number generator." },
    { MsgKey::NODT, "No time step calculation policy has been selected in the "
      "preceeding block. Use keyword 'dt' to set a constant or 'cfl' to set an "
       "adaptive time step size calculation policy." },
    { MsgKey::MULDT, "Multiple time step calculation policies has been "
      "selected in the preceeding block. Use either keyword 'dt' to set a "
      "constant or 'cfl' to set an adaptive time step size calculation policy. "
      "Setting 'cfl' and 'dt' are mutually exclusive. If both 'cfl' and 'dt' "
      "are set, 'dt' wins." },
    { MsgKey::NOINIT, "No initialization policy has been specified within the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'init' to specify an initialization policy." },
    { MsgKey::NOPROBLEM, "No test problem has been specified within the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'problem' to specify a test problem." },
    { MsgKey::NOCOEFF, "No coefficients policy has been specified within the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'coeff' to specify an coefficients policy." },
    { MsgKey::NODELTA, "No icdelta...end block with at least a single "
      "spike...end block has been specified within the block preceding this "
      "position. This is mandatory for the preceding block if the joint delta "
      "initpolicy is selected. Pick an initpolicy different than jointdelta "
      "(using keyword 'init') or specify at least a single spike...end block "
      "(within an icdelta...end block)." },
    { MsgKey::NOBETA, "No beta...end block with at least a single "
      "betapdf...end block has been specified within the block preceding this "
      "position. This is mandatory for the preceding block if jointbeta "
      "initpolicy is selected. Pick an initpolicy different than jointbeta "
      "(using keyword 'init') or specify at least a single betapdf...end block "
      "(within a icbeta...end block)." },
    { MsgKey::ODDSPIKES, "Incomplete spike...end block has been specified "
      "within the  block preceding this position. A spike...end block "
      "must contain an even number of real numbers, where every odd one is the "
      "sample space position of a spike followed by the spike height "
      "specifying the relative probability of the spike." },
    { MsgKey::WRONGBETAPDF, "Wrong number of beta distribution parameters. A "
      "beta distribution must be configured by exactly four real numbers in a "
      "betapdf...end block." },
    { MsgKey::NOTERMS, "Statistic requires at least one variable." },
    { MsgKey::NOSAMPLES, "PDF requires at least one sample space variable." },
    { MsgKey::INVALIDSAMPLESPACE, "PDF sample space specification incorrect. A "
      "non-empty list of sample space variables, must be followed by a "
      "colon, followed by a non-empty list of bin sizes (reals numbers), e.g., "
      "\"(x y : 0.1 0.2)\"" },
    { MsgKey::MALFORMEDSAMPLE, "A PDF sample space variable must be a single "
      "upper or lowercase letter optionally followed by an integer. "
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
    { MsgKey::UNFINISHED, "Block started but not finished by the 'end' "
      "keyword." },
    { MsgKey::VORTICAL_UNFINISHED, "Specifying the vortical flow test problem "
      "requires the specification of parameters alpha, beta, and p0. The error"
      "is in the block finished above the line above." },
    { MsgKey::BC_EMPTY, "Error in the preceding block. Empty boundary "
      "condition specifications, e.g., 'sideset end', are not allowed." },
    { MsgKey::WRONGSIZE, "Error in the preceding line or block. The size of "
      "the parameter vector is incorrect." },
    { MsgKey::HYDROTIMESCALES, "Error in the preceding line or block. "
      "Specification of a 'hydrotimescales' vector missing." },
    { MsgKey::HYDROPRODUCTIONS, "Error in the preceding line or block. "
      "Specification of a 'hydroproductions' vector missing." },
    { MsgKey::CHARMARG, "Arguments starting with '+' are assumed to be inteded "
      "for the Charm++ runtime system. Did you forget to prefix the command "
      "line with charmrun? If this warning persists even after running with "
      "charmrun, then Charm++ does not understand it either. See the Charm++ "
      "manual at http://charm.cs.illinois.edu/manuals/html/charm++/"
      "manual.html." }
  } );

  //! \brief Parser error and warning message handler.
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
  //! \author J. Bakosi
  template< class Stack, MsgType type, MsgKey key, class Input >
  static void Message( Stack& stack, const Input& in ) {
    const auto& msg = message.find(key);
    if (msg != message.end()) {
      std::stringstream ss;
      const std::string typestr( type == MsgType::ERROR ? "Error" : "Warning" );
      auto pos = pegtl::position_info( in );
      if (!in.empty()) {
        ss << typestr << " while parsing '" << in.string() << "' at "
           << pos.line << ',' << pos.column << ". " << msg->second;
      } else {
        ss << typestr << " while parsing at " << pos.line << ',' << pos.column
           << ". " << msg->second;
      }
      stack.template push_back< tag::error >( ss.str() );
    } else {
      stack.template push_back< tag::error >
        ( std::string("Unknown parser ") + 
          (type == MsgType::ERROR ? "error" : "warning" ) +
          " with no location information." );
    }
  }

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-local-typedef"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
  #endif
  //! \brief Compile-time test functor verifying that type U is a keyword
  //! \details This functor is used for triggering a compiler error if any of
  //!   the expected option values is not in the keywords pool of the grammar.
  //!   It is used inside of a boost::mpl::for_each to run a compile-time loop
  //!   over an MPL sequence, e.g., a vector, which verifies that each type in
  //!   the vector is a valid keyword that defines the type 'pegtl_string'.
  //!  \see kw::keyword in Control/Keyword.h
  //!  \see e.g. store_option
  // \author J. Bakosi
  template< template< class > class use >
  struct is_keyword {
    template< typename U > void operator()( U ) {
      // Attempting to define the type blow accomplishes triggering an error if
      // the type does not define pegtl_string. The compiler, however, does not
      // see that far, and generates a warning: unused type alias 'kw', so we
      // ignore it around this template.
      using kw = typename use< U >::pegtl_string;
    }
  };
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic pop
  #endif

  //! \brief Put option (i.e., a tk::Toggle) in grammar state (or stack) at a
  //!   position given by tags
  //! \details This function is used to store an option (an object deriving from
  //!   tk::Toggle) into the grammar stack. See walker::ctr::DiffEq for an
  //!   example specialization of tk::Toggle to see how an option is created
  //!   from tk::Toggle.) The grammar stack is a hiearchical tagged tuple and
  //!   the variadic list of template arguments, tags..., are used to specify
  //!   a series tags (empty structs, see Control/Tags.h) addressing a
  //!   particular field of the tagged tuple, i.e., one tag for every additional
  //!   depth level.
  //! \param[in,out] stack Grammar stack (a tagged tuple) to operate on
  //! \param[in] in Last parsed PEGTL input token
  //! \param[in] defaults Reference to a copy of the full grammar stack at the
  //!   initial state, i.e., containing the defaults for all of its fields. This
  //!   is used to detect if the user wants to overwrite an option value that
  //!   has already been set differently from the default
  //! \author J. Bakosi
  template< class Stack, template< class > class use, class Option,
            class DefaultStack, class Input, class... tags >
  static void store_option( Stack& stack,
                            const Input& in,
                            const DefaultStack& defaults ) {
    Option opt;
    auto value = in.string();
    if (opt.exist(value)) {
      auto pos = pegtl::position_info( in );
      // Emit warning on overwriting a non-default option value. This is
      // slightly inelegant. To be more elegant, we could simply call Message()
      // here, but the warning message can be more customized here (inside of
      // this function) and thus produces a more user-friendly message, compared
      // to a static message that Message() operates with due to its generic
      // nature. Instead, we emit this more user-friendly message here
      // (during parsing), instead of after parsing as part of the final parser-
      // diagnostics. We still provide location information here though.
      if (stack.template get< tags... >() != opt.value( value ) &&
          stack.template get< tags... >() != defaults.template get< tags... >())
        g_print << "\n>>> WARNING: Multiple definitions for '"
                << opt.group() << "' option. Overwriting '"
                << opt.name( stack.template get< tags... >() ) << "' with '"
                << opt.name( opt.value( value ) ) << "' at "
                << pos.line << ',' << pos.column << ".\n\n";
      stack.template set< tags... >( opt.value( value ) );
    } else {
      Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
    }
    // trigger error at compile-time if any of the expected option values
    // is not in the keywords pool of the grammar
    boost::mpl::for_each< typename Option::keywords >( is_keyword< use >() );
  }

  // Common PEGTL actions (PEGTL actions reused by multiple grammars)

  //! PEGTL action base: do nothing by default
  //! \details This base is specialized to different actions duing parsing.
  //! \author D. Frey
  template< typename Rule >
  struct action : pegtl::nothing< Rule > {};

  //! Helper for calling action::apply for multiple actions
  //! \author D. Frey
  template< typename... As >
  struct call {
    template< typename Input, typename State >
    static void apply( const Input& in, State& state ) {
        using swallow = bool[];
        (void)swallow{ ( action< As >::apply( in, state ), true )..., true };
     }
  };

  //! Rule used to trigger action(s) for a rule
  //! \author J. Bakosi
  template< class rule, class... actions >
  struct act : rule {};

  //! \author J. Bakosi
  //! \details Specialization of action for act< rule, actions... >
  template< class rule, class... actions >
  struct action< act< rule, actions... > > : call< actions... > {};

  //! Rule used to trigger action
  template< MsgType, MsgKey > struct msg : pegtl::success {};
  //! \brief Error message dispatch
  //! \details This struct and its apply function are used to dispatch a message
  //!   (e.g., error, waring) from the parser. It is simply an interface to
  //!   Message. See This struct is practically used as a functor, i.e., a
  //!   struct or class that defines the function call operator, but instead the
  //!   function call operator, PEGTL uses the apply() member function for the
  //!   actions. Thus this struct can be passed to, i.e., specialize a template,
  //!   such as tk::grm::unknown, injecting in it the desired behavior.
  //! \author J. Bakosi
  template< MsgType type, MsgKey key >
  struct action< msg< type, key > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Message< Stack, type, key >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags > struct Set : pegtl::success {};
  //! \brief Put value in state at position given by tags without conversion
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the set member function of the underlying grammar
  //!    stack, tk::Control::set.
  //! \author J. Bakosi
  template< typename tag, typename... tags >
  struct action< Set< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template set< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags > struct Store : pegtl::success {};
  //! \brief Put value in state at position given by tags with conversion
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store member function of the underlying grammar
  //!    stack, tk::Control::store.
  //! \author J. Bakosi
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
  template< typename tag, typename... tags >
  struct Store_back : pegtl::success {};
  //! Convert and push back value to vector in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack, tk::Control::store_back.
  //! \author J. Bakosi
  template< typename tag, typename...tags >
  struct action< Store_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_back< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct Store_back_back : pegtl::success {};
  //! \brief Convert and push back value to vector of back of vector in state at
  //!    position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back_back member function of the
  //!    underlying grammar stack, tk::Control::store_back_back.
  //! \author J. Bakosi
  template< typename tag, typename...tags >
  struct action< Store_back_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_back_back< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct Store_back_back_back : pegtl::success {};
  //! \brief Convert and push back value to vector of back of vector of back of
  //!   vector in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back_back_back member function of the
  //!    underlying grammar stack, tk::Control::store_back_back_back.
  //! \author J. Bakosi
  template< typename tag, typename...tags >
  struct action< Store_back_back_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_back_back_back< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename... tags >
  struct Store_switch : pegtl::success {};
  //! \brief Put true in switch in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for setting a boolean value to true in the underlying grammar
  //!    stack via the member function tk::Control::set.
  //! \author J. Bakosi
  template< typename... tags >
  struct action< Store_switch< tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template set< tags... >( true );
    }
  };

  //! Rule used to trigger action
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct store_back_option : pegtl::success {};
  //! \brief Push back option to vector in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for pushing back an option (an object deriving from
  //!   tk::Toggle) into a vector in the grammar stack. See walker::ctr::DiffEq
  //!   for an example specialization of tk::Toggle to see how an option is
  //!   created from tk::Toggle. We also do a simple sanity check here testing
  //!   if the desired option value exist for the particular option type and
  //!   error out if there is a problem. Errors and warnings are accumulated
  //!   during parsing and diagnostics are given after the parsing is finished.
  //! \author J. Bakosi
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct action< store_back_option< use, Option, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template push_back<tag,tags...>( opt.value( in.string() ) );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      boost::mpl::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct store_back_back_option : pegtl::success {};
  //! \brief Push back option to vector of back of vector in state at position
  //!   given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for pushing back an option (an object deriving from
  //!   tk::Toggle) into the back of a vector of a vector in the grammar stack.
  //!   See walker::ctr::DiffEq for an example specialization of tk::Toggle to
  //!   see how an option is created from tk::Toggle. We also do a simple sanity
  //!   check here testing if the desired option value exist for the particular
  //!   option type and error out if there is a problem. Errors and warnings are
  //!   accumulated during parsing and diagnostics are given after the parsing
  //!   is finished. This functor is similar to store_back_option but pushes the
  //!   option back to a vector of a vector.
  //! \author J. Bakosi
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct action< store_back_back_option< use, Option, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template push_back_back<tag,tags...>( opt.value( in.string() ) );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      boost::mpl::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< typename field, typename sel, typename vec,
            typename tag, typename... tags >
  struct Insert_field : pegtl::success {};
  //! \brief Convert and insert value to map at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for inserting a value into a std::map behind a key in the
  //!   underlying grammar stack via the member function
  //!   tk::Control::insert_field. We detect a recently inserted key and its
  //!   type from the companion tuple field, "selected vector", given by types,
  //!   sel and vec, and use that key to insert an associated value in a
  //!   std::map addressed by tag and tags..., requiring at least one tag to
  //!   address the map. As an example, this is used in parsing parameters
  //!   associated to a particular random number generator, such as seed.
  //!   Example input file: "mkl_mcg59 seed 2134 uniform_method accurate end".
  //!   The selected vector here is the std::vector< tk::ctr::RNGType > under
  //!   tag::sel (at the second level of the tagged tuple). The
  //!   std::vector and its member function back() are then interrogated to find
  //!   out the key type and its value (an enum value) for the particular RNG.
  //!   This key is then used to insert a new entry in the std::map under
  //!   tag::param to store the RNG parameter. Client-code is in, e.g.,
  //!   tk::rngsse::seed.
  //! \author J. Bakosi
  template< typename field, typename sel, typename vec,
            typename tag, typename...tags >
  struct action< Insert_field< field, sel, vec, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      stack.template
        insert_field< key_type, field, tag, tags... >( key, in.string() );
    }
  };

  //! Rule used to trigger action
  template< template< class > class use, class Option,
            typename field, typename sel, typename vec,
            typename tag, typename... tags >
  struct insert_option : pegtl::success {};
  //! \brief Convert and insert option value to map at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for converting and inserting an option in a std::map in the
  //!   grammar stack. An option is an object deriving from tk::Toggle. See,
  //!   e.g., walker::ctr::DiffEq for an example specialization of tk::Toggle to
  //!   see how an option is created from tk::Toggle.
  //! \author J. Bakosi
  template< template< class > class use, class Option,
            typename field, typename sel, typename vec, typename tag,
            typename... tags >
  struct action< insert_option< use, Option, field, sel, vec, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // get recently inserted key from <sel,vec>
      using key_type =
        typename Stack::template nT< sel >::template nT< vec >::value_type;
      const key_type& key = stack.template get< sel, vec >().back();
      stack.template
        insert_opt< key_type, field, typename Option::EnumType, tag, tags... >
                  ( key, Option().value(in.string()) );
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      boost::mpl::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< typename prec > struct store_precision : pegtl::success {};
  //! \brief Set numeric precision for ASCII output of floating-point values
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for setting the precision used for outputing floating-point
  //!   values into text files. We also make sure that the precision to be set
  //!   is between the correct bounds of the underlying floating-point type.
  //! \see kw::precision_info
  //! \author J. Bakosi
  template< class prec >
  struct action< store_precision< prec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using PrEx = kw::precision::info::expect;
      std::string low( in.string() );
      std::transform( begin(low), end(low), begin(low), ::tolower );
      if (low == "max") {
        const auto maxprec = PrEx::upper;
        stack.template set< tag::prec, prec >( maxprec );
      } else {
        PrEx::type precision = std::cout.precision();  // set default
        try {   //try to convert matched str to int
          precision = std::stol( in.string() );
        }
        catch ( std::exception& ) {
          Message< Stack, ERROR, MsgKey::BADPRECISION >( stack, in );
        }
        // only set precision given if it makes sense
        if (precision >= PrEx::lower && precision <= PrEx::upper)
          stack.template set< tag::prec, prec >( precision );
        else
          Message< Stack, WARNING, MsgKey::PRECISIONBOUNDS >( stack, in );
      }
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
  //! \author J. Bakosi
  template<>
  struct action< helpkw > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& cmdinfo = stack.template get< tag::cmdinfo >();
      const auto& ctrinfo = stack.template get< tag::ctrinfo >();
      auto it = cmdinfo.find( in.string() );
      if (it != cmdinfo.end()) {
        // store keyword and its info on which help was requested
        stack.template set< tag::helpkw >( { it->first, it->second, true } );
      } else {
        it = ctrinfo.find( in.string() );
        if (it != ctrinfo.end())
          // store keyword and its info on which help was requested
          stack.template set< tag::helpkw >( { it->first, it->second, false } );
        else
          Message< Stack, ERROR, MsgKey::KEYWORD >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  template< typename push > struct match_depvar : pegtl::success {};
  //! \brief Match depvar (dependent variable) to one of the selected ones
  //! \details This is used to check the set of dependent variables previously
  //!    assigned to registered differential equations (or models).
  //! \author J. Bakosi
  template< class push >
  struct action< match_depvar< push > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // convert matched string to char
      auto var = stack.template convert< char >( in.string() );
      // find matched variable in set of selected ones
      if (depvars.find( var ) != depvars.end())
        action< push >::apply( in, stack );
      else  // error out if matched var is not selected
        Message< Stack, ERROR, MsgKey::NOSUCHDEPVAR >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct match_pdfname : pegtl::success {};
  //! \brief Match PDF name to the registered ones
  //! \details This is used to check the set of PDF names dependent previously
  //!    registered to make sure all are unique.
  //! \author J. Bakosi
  template<>
  struct action< match_pdfname > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // find matched name in set of registered ones
      if (pdfnames.find( in.string() ) == pdfnames.end()) {
        pdfnames.insert( in.string() );
        stack.template
          push_back< tag::cmd, tag::io, tag::pdfnames >( in.string() );
      }
      else  // error out if name matched var is already registered
        Message< Stack, ERROR, MsgKey::PDFEXISTS >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< template< class > class use, class Option, typename sel,
            typename vec, typename... tags >
  struct check_store_option : pegtl::success {};
  //! \brief Put option in state at position given by tags if among the selected
  //! \author J. Bakosi
  template< template < class > class use, class Option, typename sel,
            typename vec, typename... tags >
  struct action< check_store_option< use, Option, sel, vec, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // error out if chosen item does not exist in selected vector
      bool exists = false;
      for (const auto& r : stack.template get< sel, vec >()) {
        if (Option().value(in.string()) == r) exists = true;
      }
      if (exists)
        action< store_back_option< use, Option, tags... > >::apply( in, stack );
      else
        Message< Stack, ERROR, MsgKey::NOTSELECTED >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct add_depvar : pegtl::success {};
  //! \brief Add depvar (dependent variable) to the selected ones
  //! \author J. Bakosi
  template<>
  struct action< add_depvar > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // convert matched string to char
      auto newvar = stack.template convert< char >( in.string() );
      // put in new dependent variable to set of already selected ones
      if (depvars.find( newvar ) == depvars.end())
        depvars.insert( newvar );
      else  // error out if depvar is already taken
        Message< Stack, ERROR, MsgKey::EXISTS >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct start_vector : pegtl::success {};
  //! \brief Start new vector in vector
  //! \author J. Bakosi
  template< class tag, class... tags >
  struct action< start_vector< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template push_back< tag, tags... >();  // no arg: use default ctor
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct start_vector_back : pegtl::success {};
  //! \brief Start new vector in back of a vector
  //! \author J. Bakosi
  template< class tag, class... tags >
  struct action< start_vector_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      // no arg: use default ctor
      stack.template push_back_back< tag, tags... >();
    }
  };

  //! Rule used to trigger action
  template< typename tk::ctr::Moment, char var = '\0' >
  struct push_term : pegtl::success {};
  //! \brief Add matched value as Term into vector of vector of statistics
  //! \author J. Bakosi
  template< tk::ctr::Moment m, char var >
  struct action< push_term< m, var > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // If var is given, push var, otherwise push first char of value
      char v(var ? var : in.string()[0]);
      // Use a shorthand of reference to vector to push_back to
      auto& stats = stack.template get< tag::stat >();
      // Push term into current vector
      stats.back().emplace_back( tk::ctr::Term( v, field, m ) );
      // If central moment, trigger mean (in statistics)
      if (m == tk::ctr::Moment::CENTRAL) {
        tk::ctr::Term term( static_cast<char>(toupper(v)),
                            field,
                            tk::ctr::Moment::ORDINARY );
        stats.insert( stats.end()-1, tk::ctr::Product( 1, term ) );
      }
      field = 0;            // reset default field
    }
  };

  //! Rule used to trigger action
  template< tk::ctr::Moment m > struct push_sample : pegtl::success {};
  //! \brief Add matched value as Term into vector of vector of PDFs
  //! \author J. Bakosi
  template< tk::ctr::Moment m >
  struct action< push_sample< m > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& pdf = stack.template get< tag::pdf >();
      // Error out if sample space already has at least 3 dimensions
      if ( pdf.back().size() >= 3 ) {
        Message< Stack, ERROR, MsgKey::MAXSAMPLES >( stack, in );
      }
      // Error out if matched sample space variable starts with a digit
      if ( std::isdigit(in.string()[0]) )
        Message< Stack, ERROR, MsgKey::MALFORMEDSAMPLE >( stack, in );
      // Push term into current vector
      pdf.back().emplace_back( tk::ctr::Term( in.string()[0], field, m ) );
      // If central moment, trigger estimation of mean (in statistics)
      if (m == tk::ctr::Moment::CENTRAL) {
        tk::ctr::Term term( static_cast<char>(toupper(in.string()[0])),
                            field,
                            tk::ctr::Moment::ORDINARY );
        auto& stats = stack.template get< tag::stat >();
        if (!stats.empty())
          stats.insert( stats.end()-1, tk::ctr::Product( 1, term ) );
        else
          stats.emplace_back( tk::ctr::Product( 1, term ) );
      }
      field = 0;            // reset default field
    }
  };

  //! Rule used to trigger action
  struct push_binsize : pegtl::success {};
  //! \brief Push matched value into vector of vector binsizes
  //! \author J. Bakosi
  template<>
  struct action< push_binsize > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& bins = stack.template get< tag::discr, tag::binsize >().back();
      // Error out if binsize vector already has at least 3 dimensions
      if ( bins.size() >= 3 ) {
        Message< Stack, ERROR, MsgKey::MAXBINSIZES >( stack, in );
      }
      // Push term into vector if larger than zero
      const auto& binsize = stack.template convert< tk::real >( in.string() );
      if ( !(binsize > std::numeric_limits< tk::real >::epsilon()) )
        Message< Stack, ERROR, MsgKey::ZEROBINSIZE >( stack, in );
      else
        bins.emplace_back( binsize );
    }
  };

  //! Rule used to trigger action
  struct push_extents : pegtl::success {};
  //! \brief Push matched value into vector of PDF extents
  //! \author J. Bakosi
  template<>
  struct action< push_extents > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& vec = stack.template get< tag::discr, tag::extent >().back();
      // Error out if extents vector already has at least 3 pairs
      if (vec.size() >= 6)
        Message< Stack, ERROR, MsgKey::MAXEXTENTS >( stack, in );
      // Error out if extents vector already has the enough pairs to match the
      // number of sample space dimensions
      if (vec.size() >=
          stack.template get< tag::discr, tag::binsize >().back().size() * 2) {
        Message< Stack, ERROR, MsgKey::INVALIDEXTENT >( stack, in );
      }
      // Push extent into vector
      vec.emplace_back( stack.template convert< tk::real >( in.string() ) );
    }
  };

  //! Rule used to trigger action
  template< class eq, class param > struct check_vector : pegtl::success {};
  //! \brief Check parameter vector
  //! \author J. Bakosi
  template< class eq, class param >
  struct action< check_vector< eq, param > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {}
  };

  //! Rule used to trigger action
  template< class eq, class param > struct check_spikes : pegtl::success {};
  //! \brief Check if the spikes parameter vector specifications are correct
  //! \details Spikes are used to specify sample-space locations and relative
  //!    probability heights for a joint-delta PDF.
  //! \author J. Bakosi
  template< class eq, class param >
  struct action< check_spikes< eq, param > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& spike =
        stack.template get< tag::param, eq, param >().back().back();
      // Error out if the number of spikes-vector is odd
      if (spike.size() % 2)
        Message< Stack, ERROR, MsgKey::ODDSPIKES >( stack, in );
      // Error out if the sum of spike heights does not add up to unity, but
      // only if the spike block is not empty (an empty spike..end block
      // is okay and is used to specify no delta spikes for a dependent
      // variable).
      if (!spike.empty()) {
        tk::real sum = 0.0;
        for (std::size_t i=1; i<spike.size(); i+=2)  // every even is a height
          sum += spike[i];
        if (std::abs(sum-1.0) > std::numeric_limits< tk::real >::epsilon())
          Message< Stack, ERROR, MsgKey::HEIGHTSPIKES >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  template< class eq, class param > struct check_betapdfs : pegtl::success {};
  //! \brief Check if the betapdf parameter vector specifications are correct
  //! \details Betapdf vectors are used to configure univariate beta
  //!   distributions.
  //! \author J. Bakosi
  template< class eq, class param >
  struct action< check_betapdfs< eq, param > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& betapdf =
        stack.template get< tag::param, eq, param >().back().back();
      // Error out if the number parameters is not four
      if (betapdf.size() != 4)
        Message< Stack, ERROR, MsgKey::WRONGBETAPDF >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct check_expectation : pegtl::success {};
  //! \brief Check if there is at least one variable in expectation
  //! \author J. Bakosi
  template<>
  struct action< check_expectation > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      if (stack.template get< tag::stat >().back().empty())
        Message< Stack, ERROR, MsgKey::NOTERMS >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct check_binsizes : pegtl::success {};
  //! \brief Check if the number of binsizes equal the PDF sample space
  //!   variables
  //! \author J. Bakosi
  template<>
  struct action< check_binsizes > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      if (stack.template get< tag::pdf >().back().size() !=
          stack.template get< tag::discr, tag::binsize >().back().size())
          Message< Stack, ERROR, MsgKey::BINSIZES >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct check_extents : pegtl::success {};
  //! \brief Check if the number of extents equal 2 * the PDF sample space
  //!    variables
  //! \author J. Bakosi
  template<>
  struct action< check_extents > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Use a shorthand to extents vector
      const auto& e = stack.template get< tag::discr, tag::extent >().back();
      // Check if the number of extents are correct
      if (!e.empty() &&
          e.size() !=
            stack.template get< tag::discr, tag::binsize >().back().size()*2)
        Message< Stack, ERROR, MsgKey::INVALIDEXTENT >( stack, in );
      // Check if the lower extents are indeed lower than the higher extents
      if (e.size() > 1 && e[0] > e[1])
        Message< Stack, ERROR, MsgKey::EXTENTLOWER >( stack, in );
      if (e.size() > 3 && e[2] > e[3])
        Message< Stack, ERROR, MsgKey::EXTENTLOWER >( stack, in );
      if (e.size() > 5 && e[4] > e[5])
        Message< Stack, ERROR, MsgKey::EXTENTLOWER >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct check_samples : pegtl::success {};
  //! \brief Check if there is at least one sample space variable in PDF
  //! \author J. Bakosi
  template<>
  struct action< check_samples > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      if (stack.template get< tag::pdf >().back().empty())
        Message< Stack, ERROR, MsgKey::NOSAMPLES >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct save_field : pegtl::success {};
  //! \brief Save field ID to parser's state so push_term can pick it up
  //! \author J. Bakosi
  template<>
  struct action< save_field > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // field ID numbers start at 0
      auto f = stack.template convert< long >( in.string() ) - 1;
      Assert( f>=0, "Field value must be non-negative in tk::grm::save_field" );
      field = static_cast< ncomp_t >( f );
    }
  };

  // Common grammar (grammar that is reused by multiple grammars)

  //! \brief Read 'token' until 'erased' trimming, i.e., not consuming, 'erased'
  //! \author J. Bakosi
  template< class token, class erased >
  struct trim :
         pegtl::seq< token,
                     pegtl::sor<
                       pegtl::until< pegtl::at< erased > >,
                       msg< ERROR, MsgKey::PREMATURE > > > {};

  //! \brief Match unknown keyword and handle error
  //! \author J. Bakosi
  template< MsgType type, MsgKey key >
  struct unknown :
         pegtl::pad< pegtl::seq< trim< pegtl::any, pegtl::space >,
                                 msg< type, key > >,
                     pegtl::blank,
                     pegtl::space > {};

  //! \brief Match alias cmdline keyword
  //! \details An alias command line keyword is prefixed by a single dash, '-'.
  //! \author J. Bakosi
  template< class keyword >
  struct alias :
         pegtl::seq<
           pegtl::one< '-' >,
           typename keyword::info::alias::type,
           pegtl::sor< pegtl::space,
                       msg< ERROR, MsgKey::ALIAS > >  > {};

  //! \brief Match verbose cmdline keyword
  //! \details A verbose command line keyword is prefixed by a double-dash,
  //!   '--'.
  //! \author J. Bakosi
  template< class keyword >
  struct verbose :
         pegtl::seq< pegtl::string<'-','-'>,
                     typename keyword::pegtl_string,
                     pegtl::space > {};

  //! \brief Read keyword 'token' padded by blank at left and space at right
  //! \author J. Bakosi
  template< class token >
  struct readkw :
         pegtl::pad< trim< token, pegtl::space >,
                     pegtl::blank,
                     pegtl::space > {};

  //! \brief Read command line 'keyword' in verbose form, i.e., '--keyword'
  //! \details This version is used if no alias is defined for the given keyword
  //! \author J. Bakosi
  template< class keyword, typename = void >
  struct readcmd :
         verbose< keyword > {};

  //! \brief Read command line 'keyword' in either verbose or alias form
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
  //! \author J. Bakosi
  template< class keyword >
  struct readcmd< keyword,
                  typename if_< false, typename keyword::info::alias >::type > :
         pegtl::sor< verbose< keyword >, alias< keyword > > {};

  //! \brief Scan input padded by blank at left and space at right and if it
  //!   matches 'keywords', apply 'actions'
  //! \details As opposed to scan_until this rule, allows multiple actions
  //! \author J. Bakosi
  template< class keyword, class... actions >
  struct scan :
           pegtl::pad< act< trim< keyword, pegtl::space >, actions... >,
                       pegtl::blank,
                       pegtl::space > {};

  //! \brief Scan input padded by blank at left and space at right and if it
  //!   matches 'keywords', apply 'action'
  //! \details This version uses an additional custom end rule. As opposed
  //!   to scan, this rule allows an additional end-rule until which parsing is
  //!   continued. The additional custom end-rule is OR'd to pegtl::space.
  //! \author J. Bakosi
  template< class keywords, class action, class end = pegtl::space >
  struct scan_until :
         pegtl::pad< act< trim< keywords, pegtl::sor< pegtl::space, end > >,
                          action >,
                     pegtl::blank,
                     pegtl::space > {};

  //! \brief Parse comment: start with '#' until eol
  //! \author J. Bakosi
  struct comment :
         pegtl::pad< trim< pegtl::one<'#'>, pegtl::eol >,
                     pegtl::blank,
                     pegtl::eol > {};

  //! \brief Ignore comments and empty lines
  //! \author J. Bakosi
  struct ignore :
         pegtl::sor< comment, pegtl::until< pegtl::eol, pegtl::space > > {};

  //! \brief Parse a number: an optional sign followed by digits
  //! \author J. Bakosi
  struct number :
         pegtl::seq< pegtl::opt< pegtl::sor< pegtl::one<'+'>,
                                             pegtl::one<'-'> > >,
                     pegtl::digit > {};

  //! \brief Plow through 'tokens' until 'endkeyword'
  //! \author J. Bakosi
  template< class endkeyword, typename... tokens >
  struct block :
         pegtl::until<
           readkw< typename endkeyword::pegtl_string >,
           pegtl::sor< comment,
                       tokens...,
                       unknown< ERROR, MsgKey::KEYWORD > > > {};

  //! \brief Plow through vector of values between keywords 'key' and
  //!   'endkeyword', calling 'insert' for each if matches and allow comments
  //!   between values
  //! \author J. Bakosi
  template< class key, class insert, class endkeyword,
            class starter, class value = number >
  struct vector :
         pegtl::seq<
           act< readkw< typename key::pegtl_string >, starter >,
           pegtl::until< readkw< typename endkeyword::pegtl_string >,
                         pegtl::sor< comment,
                                     scan< value, insert >,
                                     unknown< ERROR, MsgKey::LIST > > > > {};

  //! \brief Scan string between characters 'lbound' and 'rbound' and if matches
  //!   apply action 'insert'
  //! \author J. Bakosi
  template< class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         pegtl::if_must< pegtl::one< lbound >,
                         act< pegtl::sor< trim< pegtl::not_one< lbound >,
                                                pegtl::one< rbound > >,
                                          unknown< ERROR, MsgKey::QUOTED > >,
                              insert >,
                         pegtl::one< rbound > > {};

  //! \brief Process 'keyword' and call its 'insert' action if matches 'kw_type'
  //! \author J. Bakosi
  template< class keyword, class insert, class kw_type = pegtl::digit >
  struct process :
         pegtl::if_must<
           readkw< typename keyword::pegtl_string >,
           scan< pegtl::sor< kw_type, msg< ERROR, MsgKey::MISSING > >,
                 insert > > {};

  //! \brief Process command line 'keyword' and call its 'insert' action if
  //!   matches 'kw_type'
  //! \author J. Bakosi
  template< class keyword, class insert, class kw_type = pegtl::any >
  struct process_cmd :
         pegtl::if_must< readcmd< keyword >,
                         scan< pegtl::sor<
                                 kw_type,
                                 msg< ERROR, MsgKey::MISSING > >,
                               insert > > {};

  //! \brief Process command line switch 'keyword'
  //! \details The value of a command line switch is a boolean, i.e., it can be
  //!    either set or unset.
  //! \author J. Bakosi
  template< class keyword, typename tag, typename... tags >
  struct process_cmd_switch :
         pegtl::seq< readcmd< keyword >, Store_switch< tag, tags... > > {};

  //! \brief Generic file parser entry point: parse 'keywords' and 'ignore'
  //!   until end of file
  //! \author J. Bakosi
  template< typename keywords, typename... ign >
  struct read_file :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         ign...,
                         unknown< ERROR, MsgKey::KEYWORD > > > {};

  //! \brief Process but ignore Charm++'s charmrun arguments starting with '+'
  //! \author J. Bakosi
  struct charmarg :
         pegtl::seq< pegtl::one<'+'>,
                     unknown< WARNING, MsgKey::CHARMARG > > {};

  //! \brief Generic string parser entry point: parse 'keywords' until end of
  //!   string
  //! \author J. Bakosi
  template< typename keywords >
  struct read_string :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         charmarg,
                         unknown< ERROR, MsgKey::KEYWORD > > > {};

  //! \brief Insert RNG parameter
  //! \details A parameter here is always an option. An option is an object
  //!   deriving from tk::Toggle. See, e.g., walker::ctr::DiffEq for an example
  //!   specialization of tk::Toggle to see how an option is created from
  //!   tk::Toggle.
  //! \author J. Bakosi
  template< template< class > class use, typename keyword,
            typename option, typename field, typename sel, typename vec,
            typename... tags >
  struct rng_option :
         process< keyword,
                  insert_option< use, option, field, sel, vec, tags... >,
                  pegtl::alpha > {};

  //! \brief fieldvar: a character, denoting a variable, optionally followed by
  //!   a digit
  //! \author J. Bakosi
  template< typename var >
  struct fieldvar :
         pegtl::sor<
           pegtl::seq< var, act< pegtl::plus< pegtl::digit >, save_field > >,
           var > {};

  //! \brief term: upper or lowercase fieldvar matched to selected depvars for
  //!   stats
  //! \author J. Bakosi
  struct term :
         pegtl::sor<
           act< fieldvar< pegtl::upper >,
                match_depvar< push_term< tk::ctr::Moment::ORDINARY > > >,
           act< fieldvar< pegtl::lower >,
                match_depvar< push_term< tk::ctr::Moment::CENTRAL > > > > {};

  //! \brief sample space variable: fieldvar matched to selected depvars
  //! \author J. Bakosi
  template< class c, tk::ctr::Moment m >
  struct sample_space_var :
         scan_until<
           fieldvar< c >,
           match_depvar< push_sample< m > >,
           pegtl::one<':'> > {};

  //! \brief samples: sample space variables optionally separated by fillers
  //! \author J. Bakosi
  struct samples :
         pegtl::sor<
           sample_space_var< pegtl::upper, tk::ctr::Moment::ORDINARY >,
           sample_space_var< pegtl::lower, tk::ctr::Moment::CENTRAL >
         > {};

  //! \brief bin(sizes): real numbers as many sample space dimensions were given
  //! \author J. Bakosi
  struct bins :
         pegtl::sor<
           scan_until< number, push_binsize, pegtl::one<')'> >,
           act< pegtl::until< pegtl::at< pegtl::one<')'> >, pegtl::any >,
                msg< ERROR, MsgKey::INVALIDBINSIZE > > > {};

  //! \brief plow through expectations between characters '<' and '>'
  //! \author J. Bakosi
  struct parse_expectations :
         readkw< pegtl::seq< act< pegtl::one<'<'>, start_vector< tag::stat > >,
                             pegtl::until< pegtl::one<'>'>, term >,
                             check_expectation > > {};

  //! \brief list of sample space variables with error checking
  //! \author J. Bakosi
  struct sample_space :
         pegtl::seq<
           start_vector< tag::pdf >,
           pegtl::until< pegtl::one<':'>, samples >,
           check_samples > {};

  //! \brief extents: optional user-specified extents of PDF sample space
  //! \author J. Bakosi
  struct extents :
         pegtl::sor<
           scan_until< number, push_extents, pegtl::one<')'> >,
           act< pegtl::until< pegtl::at< pegtl::one<')'> >, pegtl::any >,
                msg< ERROR, MsgKey::INVALIDEXTENT > > > {};

  //! \brief binsizes followed by optional extents with error checking
  //! \author J. Bakosi
  struct bins_exts :
         pegtl::seq<
           start_vector< tag::discr, tag::binsize >,
           start_vector< tag::discr, tag::extent >,
           pegtl::until< pegtl::sor< pegtl::one<';'>,
                                     pegtl::at< pegtl::one<')'> > >,
                         bins >,
           pegtl::until< pegtl::one<')'>, extents >,
           check_binsizes,
           check_extents > {};

  //! \brief Match pdf description: name + sample space specification
  //! \details Example syntax (without the quotes): "name(x y z : 1.0 2.0 3.0)",
  //!    'name' is the name of the pdf, and x,y,z are sample space variables,
  //!    while 1.0 2.0 3.0 are bin sizes corresponding to the x y z sample space
  //!    dimensions, respectively.
  //! \author J. Bakosi
  struct parse_pdf :
         readkw<
           pegtl::if_must<
             act< pegtl::seq< pegtl::identifier, pegtl::at< pegtl::one<'('> > >,
                  match_pdfname >,
             pegtl::sor< pegtl::one<'('>,
                         msg< ERROR, MsgKey::KEYWORD > >,
             pegtl::sor< pegtl::seq< sample_space, bins_exts >,
                         msg< ERROR, MsgKey::INVALIDSAMPLESPACE > > > > {};

  //! \brief Match precision of floating-point numbers in digits (for text
  //!   output)
  //! \author J. Bakosi
  template< template< class > class use, class prec >
  struct precision :
         process< use< kw::precision >,
                  store_precision< prec >,
                  pegtl::alnum > {};

  //! \brief Match control parameter
  //! \author J. Bakosi
  template< typename keyword, class kw_type, typename... tags >
  struct control :
         process< keyword, Store< tags... >, kw_type > {};

  //! \brief Match discretization control parameter
  //! \author J. Bakosi
  template< typename keyword, typename Tag >
  struct discr :
         control< keyword, pegtl::digit, tag::discr, Tag > {};

  //! \brief Match component control parameter
  //! \author J. Bakosi
  template< typename keyword, typename Tag >
  struct component :
         process< keyword,
                  Store_back< tag::component, Tag >,
                  pegtl::digit > {};

  //! \brief Match interval control parameter
  //! \author J. Bakosi
  template< typename keyword, typename Tag >
  struct interval :
         control< keyword, pegtl::digit, tag::interval, Tag > {};

  //! \brief Parse statistics ... end block
  //! \author J. Bakosi
  template< template< class > class use, template< class... Ts > class store >
  struct statistics :
         pegtl::if_must< readkw< typename use< kw::statistics >::pegtl_string >,
                         block< use< kw::end >,
                                interval< use< kw::interval >,
                                          tag::stat >,
                                process< use< kw::txt_float_format >,
                                         store< tk::ctr::TxtFloatFormat,
                                                tag::flformat,
                                                tag::stat >,
                                         pegtl::alpha >,
                                precision< use, tag::stat >,
                                parse_expectations > > {};

  //! \brief Parse diagnostics ... end block
  //! \author J. Bakosi
  template< template< class > class use, template< class... Ts > class store >
  struct diagnostics :
         pegtl::if_must< readkw< typename use< kw::diagnostics >::pegtl_string >,
                         block< use< kw::end >,
                                interval< use< kw::interval >,
                                          tag::diag >,
                                process< use< kw::txt_float_format >,
                                         store< tk::ctr::TxtFloatFormat,
                                                tag::flformat,
                                                tag::diag >,
                                         pegtl::alpha >,
                                precision< use, tag::diag > > > {};

  //! \brief Match model parameter
  //! \author J. Bakosi
  template< typename keyword, typename kw_type, typename model, typename Tag >
  struct parameter :
         control< keyword, kw_type, tag::param, model, Tag > {};

  //! \brief Match rng parameter
  //! \author J. Bakosi
  template< template< class > class use, typename keyword,
            typename option, typename model, typename... tags >
  struct rng :
         process< keyword,
                  check_store_option< use,
                                      option,
                                      tag::selected,
                                      tag::rng,
                                      tag::param, model, tags... >,
                  pegtl::alpha > {};

  //! \brief Match rngs ... end block
  //! \author J. Bakosi
  template< template< class > class use, class rngs >
  struct rngblock :
         pegtl::if_must< readkw< typename use< kw::rngs >::pegtl_string >,
                         block< use< kw::end >, rngs > > {};


  //! \brief Match equation/model parameter vector
  //! \details This structure is used to match a keyword ... end block that
  //!   contains a list (i.e., a vector) of numbers. The keyword that starts the
  //!   block is passed in via the 'keyword' template argument. The 'store'
  //!   argument abstracts away a "functor" used to store the parsed values
  //!   (usually a push_back operaton on a std::vector. The 'start' argument
  //!   abstracts away the starter functor used to start the inserting operation
  //!   before parsing a value (usually a push_back on a vector using the
  //!   default value constructor). The 'check' argument abstracts away a
  //!   functor used to do error checking on the value parsed. Arguments 'eq'
  //!   and 'param' denote two levels of the hierarchy relative to tag::param,
  //!   at which the parameter vector lives. Example client-code: see
  //!   walker::deck::icbeta, or walker::deck::icdelta.
  //! \author J. Bakosi
  template< template< class > class use,
            typename keyword,
            template< class, class... > class store,
            template< class, class... > class start,
            template< class, class > class check,
            typename eq,
            typename param >
  struct parameter_vector :
         act< vector< keyword,
                      store< tag::param, eq, param >,
                      use< kw::end >,
                      start< tag::param, eq, param > >,
              check< eq, param > > {};

  //! \brief Match equation/model option vector
  //! \details This structure is used to match a keyword ... end block that
  //!   contains a list (i.e., a vector) of numbers. The keyword that starts the
  //!   block is passed in via the 'keyword' template argument. The 'store'
  //!   argument abstracts away a "functor" used to store the parsed values
  //!   (e.g. a push_back operaton on a std::vector. The 'start' argument
  //!   abstracts away the starter functor used to start the inserting operation
  //!   before parsing a value (usually a push_back on a vector using the
  //!   default value constructor). The 'check' argument abstracts away a
  //!   functor used to do error checking on the value parsed. Arguments 'eq'
  //!   and 'param' denote two levels of the hierarchy relative to tag::param,
  //!   at which the parameter vector lives. Example client-code: see
  //!   walker::deck::sde_option_vector.
  //! \author J. Bakosi
  template< template< class > class use,
            typename keyword,
            class option,
            template< class, class... > class store,
            template< class, class... > class start,
            template< class, class > class check,
            typename eq,
            typename param >
  struct option_vector :
         pegtl::if_must<
           vector<
             keyword,
             store_back_back_option< use, option, tag::param, eq, param >,
             use< kw::end >,
             start< tag::param, eq, param >,
             pegtl::alpha >,
           check< eq, param > > {};

  //! \brief Match model parameter dependent variable
  //! \author J. Bakosi
  template< template< class > class use, typename model, typename Tag >
  struct depvar :
         pegtl::if_must<
           readkw< typename use< kw::depvar >::pegtl_string >,
           scan< pegtl::sor< pegtl::alpha, msg< ERROR, MsgKey::NOTALPHA > >,
                 Store_back< tag::param, model, Tag >,
                 add_depvar > > {};

  //! \brief Match and set keyword 'title'
  //! \author J. Bakosi
  template< template< class > class use >
  struct title :
         pegtl::if_must< readkw< typename use< kw::title >::pegtl_string >,
                         quoted< Set< tag::title > > > {};

  //! \brief Match and set policy parameter
  //! \author J. Bakosi
  template< template< class > class use, typename keyword,
            typename option, typename sde, typename... tags >
  struct policy :
         process<
           keyword,
           store_back_option< use, option, tag::param, sde, tags... >,
           pegtl::alpha > {};

  //! \brief Match and set a PDF option
  //! \author J. Bakosi
  template< class keyword, class store >
  struct pdf_option :
         process< keyword, store, pegtl::alpha > {};

  //! \brief Match pdfs ... end block
  //! \author J. Bakosi
  template< template< class > class use, template< class... Ts > class store >
  struct pdfs :
         pegtl::if_must<
           tk::grm::readkw< typename use < kw::pdfs >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::interval< use< kw::interval >, tag::pdf >,
             pdf_option< use< kw::pdf_filetype >,
                         store< tk::ctr::PDFFile,
                                tag::selected,
                                tag::pdffiletype > >,
             pdf_option< use< kw::pdf_policy >,
                         store< tk::ctr::PDFPolicy,
                                tag::selected,
                                tag::pdfpolicy > >,
             pdf_option< use< kw::pdf_centering >,
                         store< tk::ctr::PDFCentering,
                                tag::selected,
                                tag::pdfctr > >,
             pdf_option< use< kw::txt_float_format >,
                         store< tk::ctr::TxtFloatFormat,
                                tag::flformat,
                                tag::pdf > >,
             precision< use, tag::pdf >,
             parse_pdf > > {};

  //! \brief Ensures that a grammar only uses keywords from a pool of
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
  //!   the keywords pool for a grammar. The pools are boost::mpl::sets and are
  //!   used to provide help on command line arguments for a given executable.
  //!   They allow compile-time iteration with boost::mpl::for_each or
  //!   generating a run-time std::map associating, e.g., keywords to their help
  //!   strings.
  //! \warning Since the default maximum number of elements in a
  //!   boost::mpl::set is 20, and increasing this number would require custom
  //!   boost::preprocessor-generated boost::mpl headers, which is cumbersome
  //!   and error-prone, instead we work with several pools here that can each
  //!   hold a maximum 20 items and use boost::mpl::or_ and boost::mpl::has_key
  //!   to search for the keyword in either of the pools. Note that
  //!   boost::mpl::or_ allows a maximum 5 template arguments by definition,
  //!   i.e., max 5 OR'd pools, which corresponds to a maximum of 5x20=100
  //!   keywords. If you need more than that, increase
  //!   BOOST_MPL_LIMIT_METAFUNCTION_ARITY at the top of this file.
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
  //! \see For example usage with a single pool, see the template typedef
  //!   walker::cmd::use in Control/Walker/CmdLine/Grammar.h and its keywords
  //!   pool, walker::ctr::CmdLine::keywords, in
  //!   Control/Walker/CmdLine/CmdLine.h.
  //! \see For example usage with multiple pools, see the template typedef,
  //!   walker::deck::use, in Control/Walker/InputDeck/Grammar.h and its
  //!   keywords pools, walker::ctr::InputDeck::keywords1,
  //!   walker::ctr::InputDeck::keywords2, etc., in
  //!   Control/Walker/InputDeck/InputDeck.h.
  //! \see http://en.cppreference.com/w/cpp/types/conditional
  //! \see http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/set.html
  //! \see http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/has_key.html
  //! \see http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/or.html
  //! \see http://www.boost.org/doc/libs/release/libs/mpl/doc/refmanual/limit-metafunction-arity.html
  //! TODO It still would be nice to generate a more developer-friendly
  //!    compiler error if the keyword is not in the pool.
  //! \author J. Bakosi
  template< typename keyword, typename pool, typename... pools >
  struct use :
         std::conditional< boost::mpl::or_<
                             boost::mpl::has_key< pool, keyword >,
                             boost::mpl::has_key< pools, keyword >... >::value,
                           keyword,
                           char >::type {};

} // grm::
} // tk::

#endif // CommonGrammar_h
