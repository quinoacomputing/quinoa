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
    EXISTS,             //!< Variable already used
    NODEPVAR,           //!< Dependent variable has not been specified
    DEPVAR_AS_MESHREF,  //!< Depvar upstream of meshref has not been specified
    LOC_NOMESHREF,      //!< Mesh location without reference mesh
    ORI_NOMESHREF,      //!< Mesh orientation without reference mesh
    MULTIMESH,          //!< If meshes are assigned, all solvers must have one
    NOSOLVE,            //!< Dependent variable to solve for has not been spec'd
    NOSUCHDEPVAR,       //!< Dependent variable has not been previously selected
    NOSUCHCOMPONENT,    //!< No such scalar component
    NOSUCHOUTVAR,       //!< Output variable label not acceptable
    NOSUCHMULTIMATVAR,  //!< Variable not acceptable for multi-material output
    POSITIVECOMPONENT,  //!< Scalar component must be positive
    NOTALPHA,           //!< Variable must be alphanumeric
    NONCOMP,            //!< No number of components selected
    LARGECOMP,          //!< Component index indexing out of max eq sys ncomp
    BADRANGE,           //!< Incorrect time range configuration
    NONMAT,             //!< No number of materials selected
    NUMMAT,             //!< Incorrect number of materials selected
    REPMATID,           //!< Repeating material id
    ONEMATID,           //!< Material id not one-based
    GAPMATID,           //!< Material id not contiguous
    NOEOS,              //!< EOS not supported
    EOSGAMMA,           //!< Wrong number of EOS gamma parameters
    EOSCV,              //!< Wrong number of EOS cv parameters
    EOSPSTIFF,          //!< Wrong number of EOS pstiff parameters
    EOSMU,              //!< Wrong number of EOS mu parameters
    EOSJWLPARAM,        //!< Wrong number of JWL EOS parameters
    EOSJWLREFSTATE,     //!< Incorrect reference state for JWL EOS
    NODT,               //!< No time-step-size policy selected
    MULDT,              //!< Multiple time-step-size policies selected
    POINTEXISTS,        //!< Point identifier already defined
    BADPRECISION,       //!< Floating point precision specification incorrect
    BOUNDS,             //!< Specified value out of bounds
    PRECISIONBOUNDS,    //!< Floating point precision spec out of bounds
    UNFINISHED,         //!< Unfinished block
    VORTICAL_UNFINISHED,//!< Vortical flow problem configuration unfinished
    ENERGY_UNFINISHED,  //!< Nonlinear energy growth problem config unfinished
    RT_UNFINISHED,      //!< Reyleigh-Taylor unstable configuration unfinished
    BC_EMPTY,           //!< Empty boundary condition block
    SYSFCTVAR,          //!< System-FCT variable index incorrect
    BGICMISSING,        //!< Background IC unspecified
    BGMATIDMISSING,     //!< Background material id unspecified
    MESHBLOCKSUPPORT,   //!< Mesh block not supported
    MESHBLOCKIDMISSING, //!< Mesh block id unspecified
    MESHBLOCKVOL,       //!< Mesh block volume unspecified
    BOXMATIDMISSING,    //!< Box material id unspecified
    BOXMATIDWRONG,      //!< Box material id incorrect
    BOXORIENTWRONG,     //!< Box orientation incorrect
    STAGBCWRONG,        //!< Stagnation BC incorrectly configured
    SKIPBCWRONG,        //!< Skip BC incorrectly configured
    SPONGEBCWRONG,      //!< Sponge BC incorrectly configured
    NONDISJOINTBC,      //!< Different BC types assigned to the same side set
    WRONGSIZE,          //!< Size of parameter vector incorrect
    WRONGMESHMOTION,    //!< Error in mesh motion dimensions
    STEADYALE,          //!< ALE + steady state not supported
    INCOMPLETEUSERFN,   //!< Incomplete user-defined function
    T0REFODD,           //!< AMR initref vector size is odd (must be even)
    T0REFNOOP,          //!< AMR t<0 refinement will be no-op
    DTREFNOOP,          //!< AMR t>0 refinement will be no-op
    PREFTOL,            //!< p-refinement tolerance out of bounds
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
    { MsgKey::EXISTS, "Dependent variable already used." },
    { MsgKey::NOSUCHDEPVAR, "Dependent variable not selected upstream in the "
      "input file. To request an output variable, a statistic, configure a PDF "
      "a involving this variable, use this variable as a coefficients policy "
      "variable, use this variable as a refinement variable, or use a "
      "dependent variable in any way, an equation must be specified upstream "
      "in the control file assigning this variable to an equation to be "
      "integrated using the depvar keyword." },
    { MsgKey::NOSUCHCOMPONENT, "Scalar component, used in conjunction with "
      "dependent variable, does not exist upstream in the input file. This "
      "happens when referring to a scalar component of a multi-component "
      "system of equations that has less than the number of total components "
      "than the one specified. Note that numbering of the components starts "
      "from 1 and their maximum value is the number specified by the 'ncomp' "
      "keyword, inclusive, if applicable for the equation block the component "
      "specification refers to. Note that there are equation system types for "
      "which the number of components are not configurable with the 'ncomp' "
      "keyword, instead their ncomp is assumed known, e.g., for compflow ncomp "
      "= 5." },
    { MsgKey::NOSUCHOUTVAR, "Scalar component label is not acceptable as a "
      "request for an output variable. Did you mean it as upper case (as a "
      "request for an instantaneous) quantity?" },
    { MsgKey::NOSUCHMULTIMATVAR, "Scalar component label is not acceptable "
      "requesting a multi-material output variable. Did you mean it as upper "
      "case (as a request for an instantaneous) quantity?" },
    { MsgKey::POSITIVECOMPONENT, "Scalar component must be positive." },
    { MsgKey::NOTALPHA, "Variable not alphanumeric." },
    { MsgKey::NODEPVAR, "Dependent variable not specified within the block "
      "preceding this position. This is mandatory for the preceding block. Use "
      "the keyword 'depvar' to specify the dependent variable." },
    { MsgKey::DEPVAR_AS_MESHREF, "Error in the preceding solver-configuration "
       "block. Dependent variable, attempted to be used as a mesh reference "
       "variable (to couple to another solver) not specified in a solver "
       "upstream. To be able to couple a solver to another one, a dependent "
       "variable of a solver, defined upstream in the input file, can be "
       "selected. This also means that the current depvar cannot be used as "
       "the mesh reference variable." },
    { MsgKey::LOC_NOMESHREF, "Location was configured without reference mesh. "
       "This is insufficient: which mesh the location should be used with? "
       "Either remove the location or add a reference mesh." },
    { MsgKey::ORI_NOMESHREF, "Orientation was configured without reference "
       "mesh. This is insufficient: which mesh the orientation should be used "
       "with? Either remove the orientation or add a reference mesh." },
    { MsgKey::MULTIMESH, "If a solver is assigned a mesh in the input/control "
       "file, all solvers must have a mesh assigned. If no solver has a mesh "
       "assigned, the (single) mesh must be specified on the command line." },
    { MsgKey::NOSOLVE, "Dependent variable to solve for not specified within "
      "the block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'solve' to specify the type of the dependent "
      "variable to solve for." },
    { MsgKey::NONCOMP, "The number of components has not been specified in the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'ncomp' to specify the number of components." },
    { MsgKey::LARGECOMP, "The component index is too large and indexes out of "
      "the total number of scalar components of the equation system "
      "configured." },
    { MsgKey::BADRANGE, "Incorrect output time range configuration. "
      "Configuration for a time range must contain exactly 3 reals, "
      "specifying mintime, maxtime, and dt, as exactly 3 reals in that order, "
      "with maxtime > mintime and 0 < dt < maxtime-mintime." },
    { MsgKey::NONMAT, "The number of materials has not been specified in the "
      "block preceding this position. This is mandatory for the preceding "
      "block. Use the keyword 'nmat' to specify the number of materials." },
    { MsgKey::NUMMAT, "The total number of materials in all the material "
      "blocks is not equal to the number of materials 'nmat' specified for "
      "this system." },
    { MsgKey::REPMATID, "Repeating material id specified in 'material ... end' "
      "block. Material ids must be unique." },
    { MsgKey::ONEMATID, "Material ids specified in 'material ... end' blocks "
      "not one-based. Material ids must begin with one." },
    { MsgKey::GAPMATID, "Material ids specified in 'material ... end' blocks "
      "have a gap. Material ids must be contiguous." },
    { MsgKey::NOEOS, "Unsupported equation of state (EOS) specified in "
      "preceding block's 'material ... end' sub-block." },
    { MsgKey::EOSGAMMA, "Incorrect number of equation of state (EOS) 'gamma' "
      "parameters configured in the preceding block's 'material ... end' "
      "sub-block. The number of components between 'gamma ... end' is "
      "incorrect, whose size must equal the number of material-ids set by "
      "keyword 'id' in that 'material ... end' sub-block." },
    { MsgKey::EOSCV, "Incorrect number of equation of state (EOS) 'cv' "
      "parameters configured in the preceding block's 'material ... end' "
      "sub-block. The number of components between 'cv... end' is "
      "incorrect, whose size must equal the number of material-ids set by "
      "keyword 'id' in that 'material ... end' sub-block." },
    { MsgKey::EOSPSTIFF, "Incorrect number of equation of state (EOS) 'pstiff' "
      "parameters configured in the preceding block's 'material ... end' "
      "sub-block. The number of components between 'pstiff ... end' "
      "is incorrect, whose size must equal the number of material-ids set by "
      "keyword 'id' in that 'material ... end' sub-block." },
    { MsgKey::EOSMU, "Incorrect number of equation of state (EOS) 'mu' "
      "parameters configured in the preceding block's 'material ... end' "
      "sub-block. The number of components between 'mu ... end' "
      "is incorrect, whose size must equal the number of material-ids set by "
      "keyword 'id' in that 'material ... end' sub-block." },
    { MsgKey::EOSJWLPARAM, "Incorrect number of JWL equation of state (EOS) "
      "parameters configured in the preceding block's 'material ... end' "
      "sub-block. The number of components between 'param ... end' "
      "is incorrect, whose size must equal the number of material-ids set by "
      "keyword 'id' in that 'material ... end' sub-block." },
    { MsgKey::EOSJWLREFSTATE, "Either reference density or reference "
      "temperature must be specified for JWL equation of state (EOS) "
      "in the preceding block's 'material ... end' sub-block. The number of "
      "components between one of these 'param ... end' is incorrect, whose "
      "size must equal the number of material-ids set by keyword 'id' in that "
      "'material ... end' sub-block." },
    { MsgKey::NODT, "No time step calculation policy has been selected in the "
      "preceeding block. Use keyword 'dt' to set a constant or 'cfl' to set an "
       "adaptive time step size calculation policy." },
    { MsgKey::MULDT, "Multiple time step calculation policies has been "
      "selected in the preceeding block. Use either keyword 'dt' to set a "
      "constant or 'cfl' to set an adaptive time step size calculation policy. "
      "Setting 'cfl' and 'dt' are mutually exclusive. If both 'cfl' and 'dt' "
      "are set, 'dt' wins." },
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
    { MsgKey::VORTICAL_UNFINISHED, "Specifying the vortical flow test problem "
      "requires the specification of parameters alpha, beta, and p0. The error"
      "is in the block finished above the line above." },
    { MsgKey::ENERGY_UNFINISHED, "Specifying the nonlinear energy growth test "
      "problem requires the specification of parameters alpha, betax, betay, "
      "betaz, ce, kappa, and r0. The error is in the block finished above the "
      "line above."},
    { MsgKey::RT_UNFINISHED, "Specifying Reyleigh-Taylor test problem "
      "requires the specification of parameters alpha, betax, betay, betaz, "
      "kappa, r0, and p0. The error is in the block finished above the line "
      "above."},
    { MsgKey::BC_EMPTY, "Error in the preceding block. Empty boundary "
      "condition specifications, e.g., 'sideset end', are not allowed." },
    { MsgKey::SYSFCTVAR, "Error in the system-FCT variable definition block. "
      "The block must list integers between 1 and 5 both inclusive." },
    { MsgKey::BGMATIDMISSING, "Error in the preceding block. "
      "The block must contain background material id." },
    { MsgKey::MESHBLOCKSUPPORT, "Error in the preceding block. "
      "Mesh block based IC not supported by discretization scheme." },
    { MsgKey::MESHBLOCKIDMISSING, "Error in the preceding block. "
      "Each IC mesh block must specify the mesh block id." },
    { MsgKey::MESHBLOCKVOL, "Error in the preceding block. "
      "Mesh block volume must be specified, if energy content is used to "
      "initialize block" },
    { MsgKey::BOXMATIDMISSING, "Error in the preceding block. "
      "Each IC box must specify material id in the box." },
    { MsgKey::BOXMATIDWRONG, "Error in the preceding block. "
      "Material id in IC box larger than number of materials." },
    { MsgKey::BOXORIENTWRONG, "Error in the preceding block. "
      "Orientation of IC box must have 3 components." },
    { MsgKey::STAGBCWRONG, "Stagnation boundary conditions incorrectly "
      "configured. Within a bc_stag ... end block there must be a point ... "
      "end block and a radius ... end block. Both point and radius blocks must "
      "contain floating-point numbers, and the number of items in the point "
      "block must be exactly 3x that of radii." },
    { MsgKey::SKIPBCWRONG, "Skip boundary conditions incorrectly "
      "configured. Within a bc_skip ... end block there must be a point ... "
      "end block and a radius ... end block. Both point and radius blocks must "
      "contain floating-point numbers, and the number of items in the point "
      "block must be exactly 3x that of radii." },
    { MsgKey::SPONGEBCWRONG, "Sponge symmetry boundary conditions incorrectly "
      "configured. Within a bc_sym ... end block, if a sponge parameter vector "
      "is given, its size must equal the number of side sets configured for "
      "symmetry BCs, and each entry must be between 0.0 and 1.0, prescribing "
      "the percentage of absorption at the boundary." },
    { MsgKey::NONDISJOINTBC, "Different boundary condition types are assigned "
      "to the same side set." },
    { MsgKey::WRONGSIZE, "Error in the preceding line or block. The size of "
      "the parameter vector is incorrect." },
    { MsgKey::WRONGMESHMOTION, "Error in the preceding line or block. Mesh "
      "motion dimension list can only involve the integers 0, 1, and 2, and "
      "the size of the list of dimensions must be lower than 4 and larger "
      "than 0." },
    { MsgKey::STEADYALE, "Error in the preceding line or block. Arbitrary "
      "Lagrangian-Eulerian mesh motion is not supported together with marching "
      "to steady state." },
    { MsgKey::INCOMPLETEUSERFN, "Error in the preceding line or block. "
      "Incomplete user-defined function. This usually means the number of "
      "entries in list is either empty (i.e., the function is not defined) or "
      "the number of entries is larger than zero but it is not divisible by "
      "the correct number. For example, if a R->R^3 function is expected the "
      "number of descrete entries must be divisible by 4: one 'column' for "
      "the abscissa, and 3 for the ordinate." },
    { MsgKey::T0REFODD, "Error in the preceding line or block. "
      "The number of edge-nodes, marking edges as pairs of nodes, used for "
      "explicit tagging of edges for initial mesh refineoment is odd (it must "
      "be even)." },
    { MsgKey::T0REFNOOP, "Initial (t<0) mesh refinement configuration will be a"
      " no-op. Initial mesh refinement requires in the amr ... end block: (1) "
      "'" +  kw::amr_t0ref::string() + " true' and at least one initial "
      "refinement type, e.g., '" + kw::amr_initial::string() + ' ' +
      kw::amr_uniform::string() + "'." },
    { MsgKey::DTREFNOOP, "Mesh refinement configuration for t>0 will be a "
      "no-op. During-timestepping (t>0) mesh refinement configuration "
      "requires in the amr ... end block: (1) '" + kw::amr_dtref::string() +
      " true' and (2) a specification of at least one refinement variable, "
      "e.g., '" + kw::amr_refvar::string() + " c end'." },
    { MsgKey::PREFTOL, "The p-refinement tolerance must be a real number "
      "between 0.0 and 1.0, both inclusive." },
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

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-local-typedef"
  #elif defined(STRICT_GNUC)
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
  #endif
  //! Compile-time test functor verifying that type U is a keyword
  //! \details This functor is used for triggering a compiler error if any of
  //!   the expected option values is not in the keywords pool of the grammar.
  //!   It is used inside of a brigand::for_each to run a compile-time loop
  //!   over an type sequence, e.g., a list, which verifies that each type in
  //!   the list is a valid keyword that defines the type 'pegtl_string'.
  //!  \see kw::keyword in Control/Keyword.h
  //!  \see e.g. store_option
  template< template< class > class use >
  struct is_keyword {
    template< typename U > void operator()( brigand::type_<U> ) {
      // Attempting to define the type below accomplishes triggering an error if
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

  #if defined(__clang__)
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wunused-template"
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
  template< class Stack, template< class > class use, class Option,
            class DefaultStack, class Input, class... tags >
  static void store_option( Stack& stack,
                            const Input& in,
                            const DefaultStack& defaults ) {
    Option opt;
    auto value = in.string();
    if (opt.exist(value)) {
      auto pos = in.position();
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
                << pos.line << ',' << pos.byte_in_line << ".\n\n";
      stack.template get< tags... >() = opt.value( value );
    } else {
      Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
    }
    // trigger error at compile-time if any of the expected option values
    // is not in the keywords pool of the grammar
    brigand::for_each< typename Option::keywords >( is_keyword< use >() );
  }
  #if defined(__clang__)
    #pragma clang diagnostic pop
  #endif

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
  template< typename tag, typename... tags > struct Set : pegtl::success {};
  //! Put value in state at position given by tags without conversion
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the set member function of the underlying grammar
  //!    stack, tk::Control::set.
  template< typename tag, typename... tags >
  struct action< Set< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >() = in.string();
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
  template< typename tag, typename... tags >
  struct Store_back : pegtl::success {};
  //! Convert and push back value to vector in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack, tk::Control::store_back.
  template< typename tag, typename...tags >
  struct action< Store_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_back< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct Store_bool : pegtl::success {};
  //! \brief Convert and push back a bool to vector of ints in state at position
  //!    given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack, tk::Control::store_back.
  template< typename tag, typename...tags >
  struct action< Store_bool< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_bool< tag, tags... >( in.string() );
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
  template< typename tag, typename...tags >
  struct action< Store_back_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template store_back_back< tag, tags... >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename tag, typename... tags >
  struct Back_back_store : pegtl::success {};
  //! \brief Convert and store value to vector of vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack. tag and tags... address a vector of vectors, whose
  //!    inner value_type is a tagged tuple to which we store here after
  //!    conversion, indexed by target.
  template< typename target, typename tag, typename...tags >
  struct action< Back_back_store< target, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().back().template
        store< target >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename tag, typename... tags >
  struct Back_store : pegtl::success {};
  //! \brief Convert and store value to vector of vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack. tag and tags... address a vector of vectors, whose
  //!    inner value_type is a tagged tuple to which we store here after
  //!    conversion, indexed by target.
  template< typename target, typename tag, typename...tags >
  struct action< Back_store< target, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().template
        store< target >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename subtarget, typename tag,
            typename... tags >
  struct Back_deep_store : pegtl::success {};
  //! \brief Convert and store value to vector of vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack. tag and tags... address a vector of vectors, whose
  //!    inner value_type is a tagged tuple to which we store here after
  //!    conversion, indexed by target and subtarget (hence deep).
  template< typename target, typename subtarget, typename tag, typename...tags >
  struct action< Back_deep_store< target, subtarget, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().template
        store< target, subtarget >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename tag, typename... tags >
  struct Back_back_store_back : pegtl::success {};
  //! \brief Convert and store value to vector of vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack. tag and tags... address a vector of vectors, whose
  //!    inner value_type is a tagged tuple to which we store_back here after
  //!    conversion, indexed by target.
  template< typename target, typename tag, typename...tags >
  struct action< Back_back_store_back< target, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().back().template
        store_back< target >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename tag, typename... tags >
  struct Back_store_back : pegtl::success {};
  //! \brief Convert and store value to vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack.
  template< typename target, typename tag, typename...tags >
  struct action< Back_store_back< target, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().template
        store_back< target >( in.string() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename subtarget, typename tag,
            typename... tags >
  struct Back_deep_store_back : pegtl::success {};
  //! \brief Convert and store value to vector of vector in state at position
  //!   given by tags and target
  //! \details This struct and its apply function are used as a functor-like
  //!    wrapper for calling the store_back member function of the underlying
  //!    grammar stack. tag and tags... address a vector of vectors, whose
  //!    inner value_type is a tagged tuple to which we store_back here after
  //!    conversion, indexed by target and subtarget (hence deep).
  template< typename target, typename subtarget, typename tag, typename...tags >
  struct action< Back_deep_store_back< target, subtarget, tag, tags... > >
  {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< tag, tags... >().back().template
        store_back< target, subtarget >( in.string() );
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
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct store_back_option : pegtl::success {};
  //! Push back option to vector in state at position given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for pushing back an option (an object deriving from
  //!   tk::Toggle) into a vector in the grammar stack. See walker::ctr::DiffEq
  //!   for an example specialization of tk::Toggle to see how an option is
  //!   created from tk::Toggle. We also do a simple sanity check here testing
  //!   if the desired option value exist for the particular option type and
  //!   error out if there is a problem. Errors and warnings are accumulated
  //!   during parsing and diagnostics are given after the parsing is finished.
  template< template < class > class use, class Option,
            typename tag, typename... tags >
  struct action< store_back_option< use, Option, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template get<tag,tags...>().push_back( opt.value( in.string() ) );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< typename target, template < class > class use,
            class Option, typename tag, typename... tags >
  struct back_back_store_option : pegtl::success {};
  //! \brief Push back option to vector of back of vector in state at position
  //!   given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for storing an option (an object deriving from tk::Toggle) in
  //!   a place in the stack. tag and tags... address a vector of vectors, whose
  //!   inner value_type is a nested tagged tuple whose field in where we store
  //!   here after conversion, indexed by target.
  //!   See walker::ctr::DiffEq for an example specialization of tk::Toggle to
  //!   see how an option is created from tk::Toggle. We also do a simple sanity
  //!   check here testing if the desired option value exist for the particular
  //!   option type and error out if there is a problem. Errors and warnings are
  //!   accumulated during parsing and diagnostics are given after the parsing
  //!   is finished.
  template< typename target, template < class > class use,
            class Option, typename tag, typename... tags >
  struct action< back_back_store_option< target, use, Option, tag, tags... > >
  {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template get< tag, tags... >().back().template
          get< target >() = opt.value( in.string() );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< typename target, template < class > class use,
            class Option, typename tag, typename... tags >
  struct back_store_option : pegtl::success {};
  //! \brief Push back option to back of vector in state at position
  //!   given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for storing an option (an object deriving from tk::Toggle) in
  //!   a place in the stack. tag and tags... address a vector, whose
  //!   inner value_type is a nested tagged tuple whose field in where we store
  //!   here after conversion, indexed by target.
  //!   See walker::ctr::DiffEq for an example specialization of tk::Toggle to
  //!   see how an option is created from tk::Toggle. We also do a simple sanity
  //!   check here testing if the desired option value exist for the particular
  //!   option type and error out if there is a problem. Errors and warnings are
  //!   accumulated during parsing and diagnostics are given after the parsing
  //!   is finished.
  template< typename target, template < class > class use,
            class Option, typename tag, typename... tags >
  struct action< back_store_option< target, use, Option, tag, tags... > >
  {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template get< tag, tags... >().back().template
          get< target >() = opt.value( in.string() );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
    }
  };

  //! Rule used to trigger action
  template< typename target, typename subtarget, template < class > class use,
            class Option, typename tag, typename... tags >
  struct back_deep_store_option : pegtl::success {};
  //! \brief Push back option to vector of back of vector in state at position
  //!   given by tags
  //! \details This struct and its apply function are used as a functor-like
  //!   wrapper for storing an option (an object deriving from tk::Toggle) in
  //!   a place in the stack. tag and tags... address a vector of vectors, whose
  //!   inner value_type is a nested tagged tuple whose field in where we store
  //!   here after conversion, indexed by target and subtarget.
  //!   See walker::ctr::DiffEq for an example specialization of tk::Toggle to
  //!   see how an option is created from tk::Toggle. We also do a simple sanity
  //!   check here testing if the desired option value exist for the particular
  //!   option type and error out if there is a problem. Errors and warnings are
  //!   accumulated during parsing and diagnostics are given after the parsing
  //!   is finished.
  template< typename target, typename subtarget, template < class > class use,
            class Option, typename tag, typename... tags >
  struct action< back_deep_store_option< target, subtarget, use, Option,
                 tag, tags... > >
  {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      Option opt;
      if (opt.exist(in.string())) {
        stack.template get< tag, tags... >().back().template
          get< target, subtarget >() = opt.value( in.string() );
      } else {
        Message< Stack, ERROR, MsgKey::NOOPTION >( stack, in );
      }
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
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
  template< template< class > class use, class Option,
            typename field, typename sel, typename vec, typename tag,
            typename... tags >
  struct action< insert_option< use, Option, field, sel, vec, tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // get recently inserted key from <sel,vec>
      const auto& key = stack.template get< sel, vec >().back();
      stack.template
        insert_field< field, typename Option::EnumType, tag, tags... >
                    ( key, Option().value(in.string()) );
      // trigger error at compile-time if any of the expected option values
      // is not in the keywords pool of the grammar
      brigand::for_each< typename Option::keywords >( is_keyword< use >() );
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
  template< class prec >
  struct action< store_precision< prec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using PrEx = kw::precision::info::expect;
      std::string low( in.string() );
      std::transform( begin(low), end(low), begin(low), ::tolower );
      if (low == "max") {
        const auto maxprec = PrEx::upper;
        stack.template get< tag::prec, prec >() = maxprec;
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
          stack.template get< tag::prec, prec >() = precision;
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
  template< typename push > struct match_depvar : pegtl::success {};
  //! \brief Match depvar (dependent variable) to one of the selected ones
  //! \details This is used to check the set of dependent variables previously
  //!    assigned to registered differential equations (or models).
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
  struct add_depvar : pegtl::success {};
  //! Add depvar (dependent variable) to the selected ones
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

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct start_vector : pegtl::success {};
  //! Start new vector in vector
  template< class tag, class... tags >
  struct action< start_vector< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template get< tag, tags... >().emplace_back();
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct start_vector_back : pegtl::success {};
  //! Start new vector in back of a vector
  template< class tag, class... tags >
  struct action< start_vector_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      // no arg: use default ctor
      stack.template get< tag, tags... >().back().emplace_back();
    }
  };

  //! Rule used to trigger action
  template< typename tag, typename... tags >
  struct start_vector_back_back : pegtl::success {};
  //! Start new vector in back of a vector
  template< class tag, class... tags >
  struct action< start_vector_back_back< tag, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      // no arg: use default ctor
      stack.template get< tag, tags... >().back().back().emplace_back();
    }
  };

  //! Rule used to trigger action
  template< class eq, class param, class... xparam >
  struct check_vector : pegtl::success {};
  //! Check parameter vector
  template< class eq, class param, class... xparam >
  struct action< check_vector< eq, param, xparam... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {}
  };

  //! Rule used to trigger action
  template< typename... tags >
  struct noop : pegtl::success {};
  //! Action that does nothing
  template< typename... tags >
  struct action< noop< tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {}
  };

  //! Rule used to trigger action
  struct save_field : pegtl::success {};
  //! Save field ID to parser's state so push_term can pick it up
  template<>
  struct action< save_field > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // field ID numbers start at 0
      auto f = stack.template convert< long >( in.string() ) - 1;
      if (f < 0)
        Message< Stack, ERROR, MsgKey::POSITIVECOMPONENT >( stack, in );
      field = static_cast< ncomp_t >( f );
    }
  };

  //! Rule used to trigger action
  template< typename Tag, typename... Tags >
  struct store_lua : pegtl::success {};
  //! Append character parsed in a lua ... end block to a string
  template< typename Tag, typename... Tags >
  struct action< store_lua< Tag, Tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      stack.template get< Tag, Tags..., tag::lua >().back() += in.string();
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

  //! \brief Scan input padded by blank at left and space at right and if it
  //!   matches 'keywords', apply 'action'
  //! \details This version uses an additional custom end rule. As opposed
  //!   to scan, this rule allows an additional end-rule until which parsing is
  //!   continued. The additional custom end-rule is OR'd to pegtl::space.
  template< class keywords, class action, class end = pegtl::space >
  struct scan_until :
         pegtl::pad< act< trim< keywords, pegtl::sor< pegtl::space, end > >,
                          action >,
                     pegtl::blank,
                     pegtl::space > {};

  //! \brief Scan input padded by blank at left and space at right and if it
  //!   exactly matches 'keyword', apply 'actions'
  template< class keyword, class... actions >
  struct exact_scan :
         pegtl::pad< act< pegtl::until< typename keyword::pegtl_string,
                                        pegtl::space >, actions... >,
                     pegtl::blank,
                     pegtl::space > {};

  //! Parse comment: start with '#' until eol
  struct comment :
         pegtl::pad< trim< pegtl::one<'#'>, pegtl::eol >,
                     pegtl::blank,
                     pegtl::eol > {};

  //! Ignore comments and empty lines
  struct ignore :
         pegtl::sor< comment, pegtl::until< pegtl::eol, pegtl::space > > {};

  //! Parse a number: an optional sign followed by digits
  struct number :
         pegtl::seq< pegtl::opt< pegtl::sor< pegtl::one<'+'>,
                                             pegtl::one<'-'> > >,
                     pegtl::digit > {};

  //! Plow through 'tokens' until 'endkeyword'
  template< class endkeyword, typename... tokens >
  struct block :
         pegtl::until<
           readkw< typename endkeyword::pegtl_string >,
           pegtl::sor< comment,
                       ignore,
                       tokens...,
                       unknown< ERROR, MsgKey::KEYWORD > > > {};

  //! \brief Read in list of dimensions between keywords 'key' and
  //!   'endkeyword', calling 'insert' for each if matches and allow comments
  //!   between values
  template< class key, class insert, class endkeyword,
            class starter = noop< key >, class value = pegtl::digit >
  struct dimensions :
         pegtl::seq<
           act< readkw< typename key::pegtl_string >, starter >,
           block< endkeyword, scan< value, insert > > > {};

  //! Plow through vector of values between keywords 'key' and
  //!   'endkeyword', calling 'insert' for each if matches and allow comments
  //!   between values
  template< class key, class insert, class endkeyword,
            class starter = noop< key >, class value = number >
  // cppcheck-suppress syntaxError
  struct vector :
         pegtl::seq<
           act< readkw< typename key::pegtl_string >, starter >,
           block< endkeyword, scan< value, insert > > > {};

  //! \brief Scan string between characters 'lbound' and 'rbound' and if matches
  //!   apply action 'insert'
  template< class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         pegtl::if_must< pegtl::one< lbound >,
                         act< pegtl::sor< trim< pegtl::not_one< lbound >,
                                                pegtl::one< rbound > >,
                                          unknown< ERROR, MsgKey::QUOTED > >,
                              insert >,
                         pegtl::one< rbound > > {};

  //! Read and store a filename between quotes
  template< template< class > class use, class tag, class... tags >
  struct filename :
          pegtl::if_must<
            readkw< typename use< kw::filename >::pegtl_string >,
            quoted< Store_back< tag, tags... > > > {};

  //! \brief Process 'keyword' and if matches, parse following token (expecting
  //!   'kw_type' and call 'insert' action on it
  template< class keyword, class insert, class kw_type = pegtl::digit >
  struct process :
         pegtl::if_must<
           readkw< typename keyword::pegtl_string >,
           scan< pegtl::sor< kw_type, msg< ERROR, MsgKey::MISSING > >,
                 insert > > {};

  //! \brief Process 'keyword' and if matches, parse following token (expecting
  //!   pegtl::alpha and call zero or more actions on it
  template< class keyword, class... actions >
  struct process_alpha :
         pegtl::if_must<
           readkw< typename keyword::pegtl_string >,
           scan< pegtl::sor< pegtl::alpha, msg< ERROR, MsgKey::MISSING > >,
                 actions... > > {};

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

  //! \brief Generic file parser entry point: parse 'keywords' and 'ignore'
  //!   until end of file
  template< typename keywords, typename... ign >
  struct read_file :
         pegtl::until< pegtl::eof,
                       pegtl::sor<
                         keywords,
                         ign...,
                         unknown< ERROR, MsgKey::KEYWORD > > > {};

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

  //! \brief Match fieldvar: a character, denoting a variable, optionally
  //!   followed by a digit
  template< typename var >
  struct fieldvar :
         pegtl::sor<
           pegtl::seq< var, act< pegtl::plus< pegtl::digit >, save_field > >,
           var > {};

  //! Match precision of floating-point numbers in digits (for text output)
  template< template< class > class use, class prec >
  struct precision :
         process< use< kw::precision >,
                  store_precision< prec >,
                  pegtl::alnum > {};

  //! Match control parameter, enforce bounds if defined
  template< typename keyword, class kw_type, template< class... > class store,
            typename... tags >
  struct control :
         pegtl::if_must<
           process< keyword, store< tags... >, kw_type >,
           typename std::conditional<
             tk::HasVar_expect_lower< typename keyword::info >::value,
             check_lower_bound< keyword, tags... >,
             pegtl::success >::type,
           typename std::conditional<
             tk::HasVar_expect_upper< typename keyword::info >::value,
             check_upper_bound< keyword, tags... >,
             pegtl::success >::type > {};

  //! Match discretization control parameter
  template< template< class > class use, typename keyword, typename Tag >
  struct discrparam :
           control< use< keyword >, pegtl::digit, Store, tag::discr, Tag > {};

  //! Match component control parameter
  template< typename keyword, typename Tag >
  struct component :
         process< keyword,
                  Store_back< tag::component, Tag >,
                  pegtl::digit > {};

  //! Match output interval control parameter in units of iteration count
  template< typename keyword, typename Tag, typename... Tags >
  struct interval_iter :
         control< keyword, pegtl::digit, Store, Tag, Tags... > {};

  //! Match output interval control parameter in units of physics time
  template< typename keyword, typename Tag, typename... Tags >
  struct interval_time :
         control< keyword, number, Store, Tag, Tags... > {};

  //! Match output range control configuration as a list of min, max, and dt
  template< template< class > class use, class keyword, class tag,
            class... tags >
  struct time_range :
         pegtl::if_must<
           tk::grm::readkw< typename use< keyword >::pegtl_string >,
           start_vector< tag, tags... >,
           tk::grm::block< use< kw::end >,
             tk::grm::scan< tk::grm::number,
               tk::grm::Store_back_back< tag, tags... > > > > {};

  //! Parse diagnostics ... end block
  template< template< class > class use, template< class... Ts > class store >
  struct diagnostics :
         pegtl::if_must< readkw< typename use< kw::diagnostics >::pegtl_string >,
                         block< use< kw::end >,
                                interval_iter< use< kw::interval_iter >,
                                  tag::output, tag::iter, tag::diag >,
                                process< use< kw::txt_float_format >,
                                         store< tk::ctr::TxtFloatFormat,
                                                tag::flformat,
                                                tag::diag >,
                                         pegtl::alpha >,
                                process< use< kw::error >,
                                         store_back_option< use,
                                                            tk::ctr::Error,
                                                            tag::diag,
                                                            tag::error >,
                                         pegtl::alpha >,
                                precision< use, tag::diag > > > {};

  //! Parse lua ... end block and store it behind Tag, Tags..., tag::lua
  template< template< class > class use, typename Tag, typename... Tags >
  struct lua :
         pegtl::if_must<
           readkw< typename use< kw::lua >::pegtl_string >,
           start_vector< Tag, Tags..., tag::lua >,
           pegtl::until< readkw< typename use< kw::end >::pegtl_string >,
                          act< pegtl::any, store_lua< Tag, Tags... > > > > {};

  //! Match model parameter
  template< typename keyword, typename kw_type, typename model, typename Tag >
  struct parameter :
         control< keyword, kw_type, Store, tag::param, model, Tag > {};

  //! Match equation/model parameter vector
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
  template< template< class > class use,
            typename keyword,
            template< class, class... > class store,
            template< class... > class start,
            template< class, class, class... > class check,
            typename eq,
            typename param,
            typename... xparams >
  struct parameter_vector :
         pegtl::if_must< vector< keyword,
                                 store< tag::param, eq, param, xparams... >,
                                 use< kw::end >,
                                 start< tag::param, eq, param, xparams... > >,
                         check< eq, param, xparams... > > {};

  //! Match model parameter dependent variable
  template< template< class > class use, typename model, typename Tag >
  struct depvar :
         pegtl::if_must<
           readkw< typename use< kw::depvar >::pegtl_string >,
           scan< pegtl::sor< pegtl::alpha, msg< ERROR, MsgKey::NOTALPHA > >,
                 Store_back< tag::param, model, Tag >,
                 add_depvar > > {};

  //! Match and set keyword 'title'
  template< template< class > class use >
  struct title :
         pegtl::if_must< readkw< typename use< kw::title >::pegtl_string >,
                         quoted< Set< tag::title > > > {};

  //! Match and set policy parameter
  template< template< class > class use, typename keyword,
            typename option, typename p, typename... tags >
  struct policy :
         process<
           keyword,
           store_back_option< use, option, tag::param, p, tags... >,
           pegtl::alpha > {};

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
