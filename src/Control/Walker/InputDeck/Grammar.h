//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 08:55:04 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's input deck grammar definition
  \details   Walker's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef WalkerInputDeckGrammar_h
#define WalkerInputDeckGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <PEGTLParsed.h>
#include <Walker/Types.h>
#include <Keywords.h>
#include <Grammar.h>

#ifdef HAS_MKL
#include <MKLGrammar.h>
#endif

#include <RNGSSEGrammar.h>

namespace walker {

extern ctr::InputDeck g_inputdeck_defaults;

//! Walker input deck facilitating user input for integrating SDEs
namespace deck {

  //! PEGTLParsed type specialized to Walker's input deck parser
  using PEGTLInputDeck =
    tk::ctr::PEGTLParsed< ctr::InputDeck,
                          pegtl::file_input< ctr::Location >,
                          tag::cmd,
                          ctr::CmdLine >;

  // Walker's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;
  //! Out-of-struct storage of field ID for pushing terms for statistics
  static int field = 0;
  //! Parser-lifetime storage for dependent variables selected. Used to track
  //! the dependent variable of differential equations assigned during parsing.
  //! It needs to be case insensitive since we only care about whether the
  //! variable is selected or not and not whether it denotes a full variable
  //! (upper case) or a fluctuation (lower case). This is true for both
  //! inserting variables into the set as well as at matching terms of products
  //! in parsing requested statistics.
  static std::set< char, tk::ctr::CaseInsensitiveCharLess > depvars;
  //! Parser-lifetime storage for PDF names. Used to track the names registere
  //! so that parsing new ones can be required to be unique.
  static std::set< std::string > pdfnames;

  // Walker's InputDeck actions

  //! start new vector in vector
  template< class tag, class... tags >
  struct start_vector : pegtl::action_base< start_vector< tag, tags... > > {
    static void apply(const std::string&, Stack& stack) {
      stack.push_back< tag, tags... >();  // no arg: use default ctor
    }
  };

  //! add matched value as Term into vector of vector of statistics
  template< tk::ctr::Moment m, char var = '\0' >
  struct push_term : pegtl::action_base< push_term< m, var > > {
    static void apply( const std::string& value, Stack& stack ) {
      // If var is given, it is triggered not user-requested
      bool plot(var ? false : true);
      // If var is given, push var, otherwise push first char of value
      char v(var ? var : value[0]);
      // Use a shorthand of reference to vector to push_back to
      auto& stats = stack.get< tag::stat >();
      // Push term into current vector
      stats.back().emplace_back( tk::ctr::Term( field, m, v, plot ) );
      // If central moment, trigger mean (in statistics)
      if (m == tk::ctr::Moment::CENTRAL) {
        tk::ctr::Term term( field, tk::ctr::Moment::ORDINARY, toupper(v), false );
        stats.insert( stats.end()-1, tk::ctr::Product( 1, term ) );
      }
      field = 0;            // reset default field
    }
  };

  //! add matched value as Term into vector of vector of PDFs
  template< tk::ctr::Moment m >
  struct push_sample : pegtl::action_base< push_sample< m > > {
    static void apply( const std::string& value, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& pdf = stack.get< tag::pdf >();
      // Error out if sample space already has at least 3 dimensions
      if ( pdf.back().size() >= 3 ) {
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::MAXSAMPLES >
                        ( stack, value );
      }
      // Error out if matched sample space variable is longer than two
      // characters or if it is two-character long but the second character is
      // not a digit - malformed sample
      if ( value.size() > 2 || (value.size() > 1 && !std::isdigit(value[1])) ) {
        tk::grm::Message< Stack,
                          tk::grm::ERROR,
                          tk::grm::MsgKey::MALFORMEDSAMPLE >
                        ( stack, value );
      }
      // Push term into current vector
      pdf.back().emplace_back( tk::ctr::Term( field, m, value[0], true ) );
      // If central moment, trigger mean (in statistics)
      if (m == tk::ctr::Moment::CENTRAL) {
        tk::ctr::Term term( field, tk::ctr::Moment::ORDINARY, toupper(value[0]),
                            false );
        auto& stats = stack.get< tag::stat >();
        if (!stats.empty())
          stats.insert( stats.end()-1, tk::ctr::Product( 1, term ) );
        else
          stats.emplace_back( tk::ctr::Product( 1, term ) );
      }
      field = 0;            // reset default field
    }
  };

  //! push matched value into vector of vector binsizes
  struct push_binsize : pegtl::action_base< push_binsize > {
    static void apply( const std::string& value, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& bins = stack.get< tag::discr, tag::binsize >().back();
      // Error out if binsize vector already has at least 3 dimensions
      if ( bins.size() >= 3 ) {
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::MAXBINSIZES >
                        ( stack, value );
      }
      // Push term into vector if larger than zero
      const auto& binsize = stack.convert< tk::real >( value );
      if ( !(binsize > std::numeric_limits< tk::real >::epsilon()) )
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::ZEROBINSIZE >
                        ( stack, value );
      else
        bins.emplace_back( binsize );
    }
  };

  //! push matched value into vector of PDF extents
  struct push_extents : pegtl::action_base< push_extents > {
    static void apply( const std::string& value, Stack& stack ) {
      // Use a shorthand of reference to vector to push_back to
      auto& vec = stack.get< tag::discr, tag::extent >().back();
      // Error out if extents vector already has at least 3 pairs
      if (vec.size() >= 6)
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::MAXEXTENTS >
                        ( stack, value );
      // Error out if extents vector already has the enough pairs to match the
      // number of sample space dimensions
      if (vec.size() >=
          stack.get< tag::discr, tag::binsize >().back().size() * 2)
        tk::grm::Message<Stack, tk::grm::ERROR, tk::grm::MsgKey::INVALIDEXTENT>
                        (stack, value);
      // Push extent into vector
      vec.emplace_back( stack.convert< tk::real >( value ) );
    }
  };

  //! check if there is at least one variable in expectation
  struct check_expectation : pegtl::action_base< check_expectation > {
    static void apply( const std::string& value, Stack& stack ) {
      if (stack.get< tag::stat >().back().empty())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOTERMS >
                        ( stack, value );
    }
  };

  //! check if the number of binsizes equal the PDF sample space variables
  struct check_binsizes : pegtl::action_base< check_binsizes > {
    static void apply( const std::string& value, Stack& stack ) {
      if (stack.get< tag::pdf >().back().size() !=
          stack.get< tag::discr, tag::binsize >().back().size())
          tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::BINSIZES >
                          ( stack, value );
    }
  };

  //! check if the number of extents equal 2 * the PDF sample space variables
  struct check_extents : pegtl::action_base< check_extents > {
    static void apply( const std::string& value, Stack& stack ) {
      // Use a shorthand to extents vector
      const auto& e = stack.get< tag::discr, tag::extent >().back();
      // Check if the number of extents are correct
      if (!e.empty() &&
          e.size() != stack.get< tag::discr, tag::binsize >().back().size()*2)
        tk::grm::Message<Stack, tk::grm::ERROR, tk::grm::MsgKey::INVALIDEXTENT>
                         (stack, value);
      // Check if the lower extents are indeed lower than the higher extents
      if (e.size() > 1 && e[0] > e[1])
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::EXTENTLOWER >
                        ( stack, value );
      if (e.size() > 3 && e[2] > e[3])
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::EXTENTLOWER >
                        ( stack, value );
      if (e.size() > 5 && e[4] > e[5])
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::EXTENTLOWER >
                        ( stack, value );
    }
  };

  //! check if there is at least one sample space variable in PDF
  struct check_samples : pegtl::action_base< check_samples > {
    static void apply( const std::string& value, Stack& stack ) {
      if (stack.get< tag::pdf >().back().empty())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOSAMPLES >
                        ( stack, value );
    }
  };

  //! save field ID to parser's state so push_term can pick it up
  struct save_field : pegtl::action_base< save_field > {
    static void apply(const std::string& value, Stack& stack) {
      field = stack.convert< int >( value ) - 1;  // field ID numbers start at 0
    }
  };

  //! add depvar (dependent variable) to the selected ones
  struct add_depvar : pegtl::action_base< add_depvar > {
    static void apply(const std::string& value, Stack& stack) {
      // convert matched string to char
      char newvar = stack.convert< char >( value );
      // put in new dependent variable to set of already selected ones
      if (depvars.find( newvar ) == depvars.end() )
        depvars.insert( newvar );
      else  // error out if depvar is already taken
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::EXISTS >
                        ( stack, value );
    }
  };

  //! match depvar (dependent variable) to one of the selected ones, used to
  //! check the set of dependent variables previously assigned to registered
  //! differential equations
  template< class push >
  struct match_depvar : pegtl::action_base< match_depvar< push > > {
    static void apply(const std::string& value, Stack& stack) {
      // convert matched string to char
      char var = stack.convert< char >( value );
      // find matched variable in set of selected ones
      if (depvars.find( var ) != depvars.end())
        push::apply( value, stack );
      else  // error out if matched var is not selected
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NODEPVAR >
                        ( stack, value );
    }
  };

  //! match PDF name to the registered ones, used to check the set of PDF names
  //! dependent previously registered to make sure all are unique
  struct match_pdfname : pegtl::action_base< match_pdfname > {
    static void apply(const std::string& value, Stack& stack) {
      // find matched name in set of registered ones
      if (pdfnames.find( value ) == pdfnames.end()) {
        pdfnames.insert( value );
        stack.push_back< tag::cmd, tag::io, tag::pdfnames >( value );
      }
      else  // error out if name matched var is already registered
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::PDFEXISTS>
                        ( stack, value );
    }
  };

  //! put option in state at position given by tags
  template< class Option, typename... tags >
  struct store_option : pegtl::action_base< store_option< Option, tags... > > {
    static void apply(const std::string& value, Stack& stack) {
      tk::grm::store_option< Stack, Option, ctr::InputDeck, tags... >
                           ( stack, value, g_inputdeck_defaults );
    }
  };

  //! put option in state at position given by tags if among the selected
  template< class Option, typename sel, typename vec, typename... tags >
  struct check_store_option :
  pegtl::action_base< check_store_option< Option, sel, vec, tags... > > {
    static void apply(const std::string& value, Stack& stack) {
      // error out if chosen item does not exist in selected vector
      bool exists = false;
      for (const auto& r : stack.template get< sel, vec >()) {
        if (Option().value(value) == r) exists = true;
      }
      if (exists) {
        tk::grm::store_back_option< Stack, Option, tags... >().apply( value,
                                                                      stack );
      } else {
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOTSELECTED >
                        ( stack, value );
      }
    }
  };

  //! set numeric precision
  struct StorePrecision : pegtl::action_base< StorePrecision > {
    static void apply(const std::string& value, Stack& stack) {
      std::string low(value);
      std::transform( begin(low), end(low), begin(low), ::tolower );
      if (low == "max") {
        stack.set< tag::discr, tag::precision >
                 ( std::numeric_limits< tk::real >::digits10 + 1 );
      } else {
        std::streamsize precision = std::cout.precision();  // set default
        try {
          precision = std::stol( value ); // try to convert matched str to int
        }
        catch ( std::exception& e ) {
          tk::grm::Message< Stack,
                            tk::grm::ERROR,
                            tk::grm::MsgKey::BADPRECISION >
                          ( stack, value );
        }
        // only set precision given if it makes sense
        if (precision > 0 &&
            precision < std::numeric_limits< tk::real >::digits10+2)
          stack.set< tag::discr, tag::precision >( precision );
        else
          tk::grm::Message< Stack,
                            tk::grm::WARNING,
                            tk::grm::MsgKey::PRECISIONBOUNDS >
                          ( stack, value );
      }
    }
  };

  // Walker's InputDeck grammar

  //! ignore: comments and empty lines
  struct ignore :
         pegtl::sor< tk::grm::comment,
                     pegtl::until< pegtl::eol, pegtl::space > > {};

  //! fieldvar: a character, denoting a variable, optionally followed by a digit
  template< typename var >
  struct fieldvar :
         pegtl::sor<
           pegtl::seq< var, pegtl::ifapply< pegtl::digit, save_field > >,
           var > {};

  //! term: upper or lowercase fieldvar matched to selected depvars for stats
  struct term :
         pegtl::sor<
           pegtl::ifapply<
             fieldvar< pegtl::upper >,
             match_depvar< push_term< tk::ctr::Moment::ORDINARY > > >,
           pegtl::ifapply<
             fieldvar< pegtl::lower >,
             match_depvar< push_term< tk::ctr::Moment::CENTRAL > > > > {};

  //! sample space variable: fieldvar matched to selected depvars
  template< class c, tk::ctr::Moment m >
  struct sample_space_var :
         tk::grm::scan_until< fieldvar< c >,
                              match_depvar< push_sample< m > >,
                              pegtl::one<':'> > {};

  //! sample space variables optionally separated by fillers
  struct samples :
         pegtl::sor< sample_space_var< pegtl::upper, tk::ctr::Moment::ORDINARY >,
                     sample_space_var< pegtl::lower, tk::ctr::Moment::CENTRAL > >
         {};

  //! bin(sizes): real numbers as many sample space dimensions were given
  struct bins :
         pegtl::sor< tk::grm::scan_until< tk::grm::number,
                                          push_binsize,
                                          pegtl::one<')'> >,
                     pegtl::ifapply<
                       pegtl::until< pegtl::at< pegtl::one<')'> >, pegtl::any >,
                       tk::grm::error< Stack,
                                       tk::grm::MsgKey::INVALIDBINSIZE > > > {};

  //! plow through expectations between characters '<' and '>'
  struct parse_expectations :
         tk::grm::readkw<
           pegtl::ifmust< pegtl::one<'<'>,
                          pegtl::apply< start_vector< tag::stat > >,
                          pegtl::until< pegtl::one<'>'>, term >,
                          pegtl::apply< check_expectation > > > {};

  //! list of sample space variables with error checking
  struct sample_space :
         pegtl::seq<
           pegtl::apply< start_vector< tag::pdf > >,
           pegtl::until< pegtl::one<':'>, samples >,
           pegtl::apply< check_samples > > {};

  //! extents: optional user-specified extents of PDF sample space
  struct extents :
         pegtl::sor< tk::grm::scan_until< tk::grm::number,
                                          push_extents,
                                          pegtl::one<')'> >,
                     pegtl::ifapply<
                       pegtl::until< pegtl::at< pegtl::one<')'> >, pegtl::any >,
                       tk::grm::error< Stack,
                                       tk::grm::MsgKey::INVALIDEXTENT > > > {};

  //! binsizes followed by optional extents with error checking
  struct bins_exts :
         pegtl::seq<
           pegtl::apply< start_vector< tag::discr, tag::binsize >,
                         start_vector< tag::discr, tag::extent > >,
           pegtl::until< pegtl::sor< pegtl::one<';'>,
                                     pegtl::at< pegtl::one<')'> > >, bins >,
           pegtl::until< pegtl::one<')'>, extents >,
           pegtl::apply< check_binsizes >,
           pegtl::apply< check_extents > > {};

  //! PDF name: a C-language identifier (alphas, digits and underscores, no
  //! leading digit), matched to already selected pdf name requiring unique
  //! names
  struct pdf_name :
         pegtl::ifapply< pegtl::identifier, match_pdfname > {};

  //! match sample space specification between characters '(' and ')', example
  //! syntax (without the quotes): "(x y z : 1.0 2.0 3.0)", where x,y,z are
  //! sample space variables, while 1.0 2.0 3.0 are bin sizes corresponding to
  //! the x y z sample space dimensions, respectively
  struct parse_pdf :
         tk::grm::readkw<
           pegtl::ifmust<
             pegtl::seq< pdf_name, pegtl::at< pegtl::one<'('> > >,
             pegtl::sor< pegtl::one<'('>,
                         pegtl::apply<
                           tk::grm::error< Stack,
                                           tk::grm::MsgKey::KEYWORD > > >,
             pegtl::sor< pegtl::seq< sample_space, bins_exts >,
                         pegtl::apply<
                           tk::grm::error< Stack,
                           tk::grm::MsgKey::INVALIDSAMPLESPACE > > > > > {};

  //! precision of floating-point numbers in digits (for text output)
  struct precision :
         tk::grm::process< Stack,
                           kw::precision::pegtl_string,
                           StorePrecision,
                           pegtl::alnum > {};

  //! control parameter
  template< typename keyword, class kw_type, typename... tags >
  struct control :
         tk::grm::process< Stack,
                           typename keyword::pegtl_string,
                           tk::grm::Store< Stack, tags... >,
                           kw_type > {};

  //! discretization control parameter
  template< typename keyword, typename Tag >
  struct discr :
         control< keyword, pegtl::digit, tag::discr, Tag > {};

  //! component control parameter
  template< typename keyword, typename Tag >
  struct component :
         tk::grm::process< Stack,
                           typename keyword::pegtl_string,
                           tk::grm::Store_back< Stack, tag::component, Tag >,
                           pegtl::digit > {};

  //! interval control parameter
  template< typename keyword, typename Tag >
  struct interval :
         control< keyword, pegtl::digit, tag::interval, Tag > {};

  //! model parameter
  template< typename keyword, typename kw_type, typename model, typename Tag >
  struct parameter :
         control< keyword, kw_type, tag::param, model, Tag > {};

  //! model parameter vector
  template< typename keyword, typename...tags >
  struct parameter_vector :
         tk::grm::vector<
           Stack,
           typename keyword::pegtl_string,
           tk::grm::Store_back_back< Stack, tag::param, tags... >,
           kw::end,
           pegtl::apply< start_vector< tag::param, tags... > > > {};

  //! rng parameter
  template< typename keyword, typename option, typename model,
            typename... tags >
  struct rng :
         tk::grm::process< Stack,
                           typename keyword::pegtl_string,
                           check_store_option< option,
                                               tag::selected,
                                               tag::rng,
                                               tag::param, model, tags... >,
                           pegtl::alpha > {};

  //! scan selected option
  template< typename keyword, typename option, typename... tags >
  struct select_option :
         tk::grm::scan< typename keyword::pegtl_string,
                        store_option< option, tag::selected, tags... > > {};

  //! scan selected option and trigger
  template< typename keyword, typename option, typename Tag,
            typename... triggers >
  struct select_option_and_trigger :
         tk::grm::scan< typename keyword::pegtl_string,
                        store_option< option, tag::selected, Tag >,
                        triggers... > {};

  //! model parameter depvar (dependent variable)
  template< typename model, typename Tag >
  struct depvar :
         pegtl::ifmust<
           tk::grm::readkw< kw::depvar::pegtl_string >,
           tk::grm::scan<
             pegtl::sor< pegtl::alpha,
                         pegtl::apply<
                           tk::grm::error< Stack, tk::grm::MsgKey::NOTALPHA > > >,
             tk::grm::Store_back< Stack, tag::param, model, Tag >,
             add_depvar > > {};

  //! scan and store_back sde keyword and option
  template< typename keyword >
  struct scan_sde :
         tk::grm::scan< typename keyword::pegtl_string,
                        tk::grm::store_back_option< Stack,
                                                    ctr::DiffEq,
                                                    tag::selected,
                                                    tag::diffeq > > {};

  //! title
  struct title :
         pegtl::ifmust< tk::grm::readkw< kw::title::pegtl_string >,
                                         tk::grm::quoted<
                                           Stack,
                                           tk::grm::Set< Stack,
                                                         tag::title > > > {};

  //! PDF option
  template< class keyword, class Option, class Tag >
  struct pdf_option :
   tk::grm::process< Stack,
                     typename keyword::pegtl_string,
                     store_option< Option, tag::selected, Tag >,
                     pegtl::alpha > {};

  //! statistics block
  struct statistics :
         pegtl::ifmust< tk::grm::readkw< kw::statistics::pegtl_string >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        interval< kw::interval, tag::stat >,
                                        parse_expectations > > {};
  //! pdfs block
  struct pdfs :
         pegtl::ifmust<
           tk::grm::readkw< kw::pdfs::pegtl_string >,
           tk::grm::block<
             Stack,
             kw::end,
             interval< kw::interval, tag::pdf >,
             pdf_option< kw::pdf_filetype, tk::ctr::PDFFile, tag::pdffiletype >,
             pdf_option< kw::pdf_policy, tk::ctr::PDFPolicy, tag::pdfpolicy >,
             pdf_option< kw::pdf_centering, tk::ctr::PDFCentering, tag::pdfctr >,
             pdf_option< kw::txt_float_format,
                         tk::ctr::TxtFloatFormat,
                         tag::float_format >,
             precision,
             parse_pdf > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< discr< kw::npar, tag::npar >,
                     discr< kw::nstep, tag::nstep >,
                     discr< kw::term, tag::term >,
                     discr< kw::dt, tag::dt >,
                     interval< kw::ttyi, tag::tty > > {};

  //! rngs
  struct rngs :
         pegtl::sor<
                     #ifdef HAS_MKL
                     tk::mkl::rngs< Stack,
                                    tag::selected, tag::rng,
                                    tag::param, tag::rngmkl >,
                     #endif
                     tk::rngsse::rngs< Stack,
                                       tag::selected, tag::rng,
                                       tag::param, tag::rngsse > > {};

  // RNGs block
  struct rngblock :
         pegtl::ifmust< tk::grm::readkw< kw::rngs::pegtl_string >,
                        tk::grm::block< Stack, kw::end, rngs > > {};

  //! policy parameter
  template< typename keyword, typename option, typename sde, typename... tags >
  struct policy :
         tk::grm::process<
           Stack,
           typename keyword::pegtl_string,
           tk::grm::store_back_option< Stack, option, tag::param, sde, tags... >,
           pegtl::alpha > {};


  //! Diagonal Ornstein-Uhlenbeck SDE
  struct diag_ornstein_uhlenbeck :
         pegtl::ifmust< scan_sde< kw::diag_ornstein_uhlenbeck >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::diagou, tag::depvar >,
                                        component< kw::ncomp, tag::diagou >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::diagou,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::diagou,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::diagou,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_sigma,
                                                          tag::diagou,
                                                          tag::sigma >,
                                        parameter_vector< kw::sde_theta,
                                                          tag::diagou,
                                                          tag::theta >,
                                        parameter_vector< kw::sde_mu,
                                                          tag::diagou,
                                                          tag::mu > > > {};

  //! Ornstein-Uhlenbeck SDE
  struct ornstein_uhlenbeck :
         pegtl::ifmust< scan_sde< kw::ornstein_uhlenbeck >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::ou, tag::depvar >,
                                        component< kw::ncomp, tag::ou >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::ou,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::ou,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::ou,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_sigma,
                                                          tag::ou,
                                                          tag::sigma >,
                                        parameter_vector< kw::sde_theta,
                                                          tag::ou,
                                                          tag::theta >,
                                        parameter_vector< kw::sde_mu,
                                                          tag::ou,
                                                          tag::mu > > > {};

  //! Skew-normal SDE
  struct skewnormal :
         pegtl::ifmust< scan_sde< kw::skewnormal >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::skewnormal, tag::depvar >,
                                        component< kw::ncomp, tag::skewnormal >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::skewnormal,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::skewnormal,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::skewnormal,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_T,
                                                          tag::skewnormal,
                                                          tag::timescale >,
                                        parameter_vector< kw::sde_sigma,
                                                          tag::skewnormal,
                                                          tag::sigma >,
                                        parameter_vector< kw::sde_lambda,
                                                          tag::skewnormal,
                                                          tag::lambda > > > {};

  //! Beta SDE
  struct beta :
         pegtl::ifmust< scan_sde< kw::beta >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::beta, tag::depvar >,
                                        component< kw::ncomp, tag::beta >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::beta,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::beta,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::beta,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_b,
                                                          tag::beta,
                                                          tag::b >,
                                        parameter_vector< kw::sde_S,
                                                          tag::beta,
                                                          tag::S >,
                                        parameter_vector< kw::sde_kappa,
                                                          tag::beta,
                                                          tag::kappa > > > {};

  //! Gamma SDE
  struct gamma :
         pegtl::ifmust< scan_sde< kw::gamma >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::gamma, tag::depvar >,
                                        component< kw::ncomp, tag::gamma >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::gamma,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::gamma,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::gamma,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_b,
                                                          tag::gamma,
                                                          tag::b >,
                                        parameter_vector< kw::sde_S,
                                                          tag::gamma,
                                                          tag::S >,
                                        parameter_vector< kw::sde_kappa,
                                                          tag::gamma,
                                                          tag::kappa > > > {};

  //! Dirichlet SDE
  struct dirichlet :
         pegtl::ifmust< scan_sde< kw::dirichlet >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::dirichlet, tag::depvar >,
                                        component< kw::ncomp, tag::dirichlet >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::dirichlet,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::dirichlet,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::dirichlet,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_b,
                                                          tag::dirichlet,
                                                          tag::b >,
                                        parameter_vector< kw::sde_S,
                                                          tag::dirichlet,
                                                          tag::S >,
                                        parameter_vector< kw::sde_kappa,
                                                          tag::dirichlet,
                                                          tag::kappa > > > {};
  //! Generalized Dirichlet SDE
  struct generalized_dirichlet :
         pegtl::ifmust< scan_sde< kw::gendir >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::gendir, tag::depvar >,
                                        component< kw::ncomp, tag::gendir >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::gendir,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::gendir,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::gendir,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_b,
                                                          tag::gendir,
                                                          tag::b >,
                                        parameter_vector< kw::sde_S,
                                                          tag::gendir,
                                                          tag::S >,
                                        parameter_vector< kw::sde_kappa,
                                                          tag::gendir,
                                                          tag::kappa >,
                                        parameter_vector< kw::sde_c,
                                                          tag::gendir,
                                                          tag::c > > > {};

  //! Wright-Fisher SDE
  struct wright_fisher :
         pegtl::ifmust< scan_sde< kw::wrightfisher >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        depvar< tag::wrightfisher, tag::depvar >,
                                        component< kw::ncomp, tag::wrightfisher >,
                                        rng< kw::rng,
                                             tk::ctr::RNG,
                                             tag::wrightfisher,
                                             tag::rng >,
                                        policy< kw::init,
                                                tk::ctr::InitPolicy,
                                                tag::wrightfisher,
                                                tag::initpolicy >,
                                        policy< kw::coeff,
                                                tk::ctr::CoeffPolicy,
                                                tag::wrightfisher,
                                                tag::coeffpolicy >,
                                        parameter_vector< kw::sde_omega,
                                                          tag::wrightfisher,
                                                          tag::omega > > > {};

  //! stochastic differential equations
  struct sde :
         pegtl::sor< dirichlet,
                     generalized_dirichlet,
                     wright_fisher,
                     ornstein_uhlenbeck,
                     diag_ornstein_uhlenbeck,
                     skewnormal,
                     gamma,
                     beta > {};

  //! 'walker' block
  struct walker :
         pegtl::ifmust< tk::grm::readkw< kw::walker::pegtl_string >,
                        tk::grm::block< Stack,
                                        kw::end,
                                        discretization_parameters,
                                        sde,
                                        rngblock,
                                        statistics,
                                        pdfs > > {};

  //! main keywords
  struct keywords :
         pegtl::sor< title, walker > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< Stack, keywords, ignore > {};

} // deck::
} // walker::

#endif // WalkerInputDeckGrammar_h
