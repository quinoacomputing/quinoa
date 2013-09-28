//******************************************************************************
/*!
  \file      src/Control/QuinoaGrammar.h
  \author    J. Bakosi
  \date      Fri Sep 27 18:16:02 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa grammar definition
  \details   Grammar definition. We use the Parsing Expression Grammar Template
             Library (PEGTL) to create the grammar and the associated parser.
             Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word of
             advice: read from the bottom up.
*/
//******************************************************************************
#ifndef QuinoaGrammar_h
#define QuinoaGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <QuinoaControlTypes.h>
#include <Option.h>
#include <QuinoaKeywords.h>

namespace quinoa {
//! Grammar definition: state, actions, grammar
namespace grm {

  using namespace pegtl;

  // State

  //! Everything is stored in Stack during parsing
  using Stack = QuinoaControl;
  //! Out-of-struct storage of field ID for pushing terms for statistics
  static int field = 0;

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

  // Actions

  //! error handler
  template< Error key >
  struct error : action_base< error<key> > {
    static void apply(const std::string& value, Stack& stack) {
      const auto& msg = err_msg.find(key);
      Assert(msg != err_msg.end(), ExceptType::FATAL,
            "Unknown parser error. You are doing something pretty nasty, e.g., "
            "overriding compiler errors by casting some arbitrary type into an "
            "enum class. Please don't rape the language. There is a good "
            "reason we rely on strong type checking.");
      if (msg != err_msg.end()) {
        Throw(ExceptType::FATAL,
              "Error while parsing '" + value + "'. " + msg->second + ".");
      }
      IGNORE(stack);
    }
  };

  //! put value in state at position given by tags without conversion
  template< typename... tags >
  struct set : action_base< set<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.set<tags...>(value);
    }
  };

  //! convert and put value in state at position given by tags
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store<tags...>(value);
    }
  };

  //! convert and push back value to vector in state at position given by tags
  template< typename...tags >
  struct store_back : action_base< store_back<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store_back<tags...>(value);
    }
  };

  //! start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value, Stack& stack) {
      stack.push_back<ctr::stats>(ctr::Product());
      IGNORE(value);        // suppress compiler warning on unused variable
    }
  };

  //! add matched value as Term into vector of Product in vector of statistics
  template< ctr::Quantity q, ctr::Moment m, char name='\0' >
  struct push_term : action_base< push_term<q, m, name> > {
    static void apply(const std::string& value, Stack& stack) {
      // If name is given, push name, otherwise push first char of value
      char na(name ? name : value[0]);
      // If name is given, it is triggered not user-requested
      bool plot(name ? false : true);
      // Use stats for shorthand of reference to stats vector
      std::vector<ctr::Product>& stats = stack.get<ctr::stats>();
      // Push term into current product
      stats.back().push_back(ctr::Term(field, q, m, na, plot));
      // If central moment, trigger mean
      if (m == ctr::Moment::CENTRAL) {
        ctr::Term term(field, q, ctr::Moment::ORDINARY, toupper(na), false);
        stats.insert(stats.end()-1, ctr::Product(1,term));
      }
      field = 0;            // reset default field
    }
  };

  //! save field ID so push_term can pick it up
  struct save_field : action_base< save_field > {
    static void apply(const std::string& value, Stack& stack) {
      field = stack.convert<int>(value) - 1;  // field ID numbering start from 0
    }
  };

  //! convert and put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : action_base< store_option<OptionType, tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      ctr::Option<OptionType> opt;
      //! Emit warning on overwrite
      if (stack.get<tags...>() != QuinoaDefaults.get<tags...>()) {
        std::cout << "\n>>> PARSER WARNING: Multiple definitions for '"
                  << opt.group() << "' option. Overwriting '"
                  << opt.name(stack.get<tags...>()) << "' with '"
                  << opt.name(opt.value(value)) << "'.\n\n";
      }
      stack.set<tags...>(opt.value(value));
    }
  };

  // Grammar

  //! read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  //! read 'token' padded by blank at left and space at right
  template< class token >
  struct read :
         pad< trim<token, space>, blank, space > {};

  //! parse input padded by blank at left and space at right and if it matches
  //! 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions... >, blank, space > {};

  // match unknown keyword and handle error
  template< Error key >
  struct unknown :
         pad< ifapply< trim<any, space>, error<key> >, blank, space > {};

  //! comment: start with '#' until eol
  struct comment :
         pad< trim<one<'#'>,eol>, blank, eol> {};

  //! plow through 'tokens' until 'endkeyword'
  template< typename endkeyword, typename... tokens >
  struct block :
         until< read<endkeyword>,
                sor<comment, tokens..., unknown<Error::KEYWORD>> > {};

  //! rng: one of the random number generators
  struct rng :
         sor< kw::mkl_mcg31::pegtl_string,
              kw::mkl_r250::pegtl_string,
              kw::mkl_mrg32k3a::pegtl_string,
              kw::mkl_mcg59::pegtl_string,
              kw::mkl_wh::pegtl_string,
              kw::mkl_mt19937::pegtl_string,
              kw::mkl_mt2203::pegtl_string,
              kw::mkl_sfmt19937::pegtl_string,
              kw::mkl_sobol::pegtl_string,
              kw::mkl_niederr::pegtl_string,
              kw::mkl_iabstract::pegtl_string,
              kw::mkl_dabstract::pegtl_string,
              kw::mkl_sabstract::pegtl_string,
              kw::mkl_nondeterm::pegtl_string > {};

  //! number: optional sign followed by digits
  struct number :
         seq< opt< sor<one<'+'>, one<'-'>> >, digit> {};

  //! plow through list of values between keywords 'key' and "end", calling
  //! 'insert' for each if matches and allowing comments between values
  template< class key, class insert, class value = number >
  struct list :
         ifmust< read<key>,
                 until< read<kw::end::pegtl_string>,
                        sor<comment,
                            parse<value,insert>,
                            unknown<Error::LIST>> > > {};

  //! scan string between characters 'lbound' and 'rbound' and if matches apply
  //! action 'insert'
  template< class insert, char lbound = '"', char rbound = '"' >
  struct quoted :
         ifmust< one<lbound>,
                 ifapply< sor<trim<not_one<lbound>, one<rbound>>,
                              unknown<Error::QUOTED>>, insert >,
                 one<rbound>
               > {};

  //! process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<keywords,insert> > {};

  //! process 'keyword' and call its 'insert' action for string matched
  //! between characters 'lbound' and 'rbound'
  template< class keyword, class insert, char lbound='"', char rbound='"' >
  struct process_quoted :
         ifmust< read<keyword>,
                 sor< quoted<insert,lbound,rbound>, unknown<Error::QUOTED>> > {};

  //! moment: 'keyword' optionally followed by a digit, pushed to vector of terms
  template< class keyword, ctr::Quantity q, ctr::Moment m >
  struct moment :
         sor < ifapply<seq<keyword, ifapply<digit,save_field>>, push_term<q,m>>,
               ifapply<keyword, push_term<q,m>> > {};

  //! terms recognized within an expectation and their mapping
  struct terms :
         sor< moment<kw::transported_scalar::pegtl_string,
                     ctr::Quantity::SCALAR,
                     ctr::Moment::ORDINARY>,
              moment<kw::transported_scalar_fluctuation::pegtl_string,
                     ctr::Quantity::SCALAR,
                     ctr::Moment::CENTRAL>,
              moment<kw::velocity_x::pegtl_string,
                     ctr::Quantity::VELOCITY_X,
                     ctr::Moment::ORDINARY>,
              unknown<Error::MOMENT>
            > {};

  //! plow through terms in expectation until character 'rbound'
  template< char rbound >
  struct expectation :
         until< one<rbound>, terms > {};

  //! plow through expectations between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_expectations :
         read< ifmust< one<lbound>, apply<start_product>, expectation<rbound> >
             > {};

  //! title
  struct title :
         ifmust< read<kw::title::pegtl_string>, quoted<set<ctr::title>> > {};

  //! analytic_geometry block
  struct analytic_geometry:
         ifmust< parse< kw::analytic_geometry::pegtl_string,
                        store_option<sel::Geometry,
                                     ctr::selected,
                                     ctr::geometry> >,
                 block<kw::end::pegtl_string> > {};

  //! discrete_geometry block
  struct discrete_geometry:
         ifmust< parse< kw::discrete_geometry::pegtl_string,
                        store_option<sel::Geometry,
                                     ctr::selected,
                                     ctr::geometry> >,
                 block<kw::end::pegtl_string,
                       process_quoted<kw::input::pegtl_string,
                                      set<ctr::io,ctr::input>>,
                       process_quoted<kw::output::pegtl_string,
                                      set<ctr::io, ctr::geomoutput>> > > {};

  //! dir block
  struct dir :
         ifmust< parse< kw::mix_dir::pegtl_string,
                        store_option<sel::Mix, ctr::selected, ctr::mix> >,
                 block< kw::end::pegtl_string,
                        process<kw::nscalar::pegtl_string,
                                store<ctr::component, ctr::nscalar>>,
                        list< kw::dir_B::pegtl_string,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::b> >,
                        list< kw::dir_S::pegtl_string,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::S> >,
                        list< kw::dir_kappa::pegtl_string,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::kappa> > > > {};

  //! gendir block
  struct gendir :
         ifmust< parse< kw::mix_gendir::pegtl_string,
                        store_option<sel::Mix, ctr::selected, ctr::mix> >,
                 block< kw::end::pegtl_string,
                        process< kw::nscalar::pegtl_string,
                                 store<ctr::component, ctr::nscalar> >,
                        list< kw::dir_B::pegtl_string,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::b> >,
                        list< kw::dir_S::pegtl_string,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::S> >,
                        list< kw::dir_kappa::pegtl_string,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::kappa> >,
                        list< kw::gendir_C::pegtl_string,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::c>> > > {};

  //! statistics block
  struct statistics :
         ifmust< read<kw::statistics::pegtl_string>,
                 block<kw::end::pegtl_string, parse_expectations<'<','>'>> > {};

  //! slm block
  struct slm :
         ifmust< parse< kw::hydro_slm::pegtl_string,
                        store_option<sel::Hydro, ctr::selected, ctr::hydro>,
                        // trigger estimating the diagonal of Reynolds-stress
                        start_product,
                        push_term<ctr::Quantity::VELOCITY_X,
                                  ctr::Moment::CENTRAL, 'u'>,
                        push_term<ctr::Quantity::VELOCITY_X,
                                  ctr::Moment::CENTRAL, 'u'>,
                        start_product,
                        push_term<ctr::Quantity::VELOCITY_Y,
                                  ctr::Moment::CENTRAL, 'v'>,
                        push_term<ctr::Quantity::VELOCITY_Y,
                                  ctr::Moment::CENTRAL, 'v'>,
                        start_product,
                        push_term<ctr::Quantity::VELOCITY_Z,
                                  ctr::Moment::CENTRAL, 'w'>,
                        push_term<ctr::Quantity::VELOCITY_Z,
                                  ctr::Moment::CENTRAL, 'w'>>,
                 block< kw::end::pegtl_string,
                        process<kw::SLM_C0::pegtl_string,
                                store<ctr::param, ctr::slm, ctr::c0>>,
                        process<kw::nvelocity::pegtl_string,
                                store<ctr::component, ctr::nvelocity>>
                      > > {};

  //! freq_gamma block
  struct freq_gamma :
         ifmust< parse< kw::freq_gamma::pegtl_string,
                        store_option<sel::Frequency,
                                     ctr::selected,
                                     ctr::frequency> >,
                 block< kw::end::pegtl_string,
                        process<kw::nfreq::pegtl_string,
                                store<ctr::component, ctr::nfrequency>>,
                        process<kw::freq_gamma_C1::pegtl_string,
                                store<ctr::param, ctr::gamma, ctr::c1>>,
                        process<kw::freq_gamma_C2::pegtl_string,
                                store<ctr::param, ctr::gamma, ctr::c2>>,
                        process<kw::freq_gamma_C3::pegtl_string,
                                store<ctr::param, ctr::gamma, ctr::c3>>,
                        process<kw::freq_gamma_C4::pegtl_string,
                                store<ctr::param, ctr::gamma, ctr::c4>> >
               > {};

  //! beta block
  struct beta :
         ifmust< parse<kw::mass_beta::pegtl_string,
                       store_option<sel::Mass, ctr::selected, ctr::mass>>,
                 block< kw::end::pegtl_string,
                        process<kw::ndensity::pegtl_string,
                                store<ctr::component, ctr::ndensity>>,
                        process<kw::Beta_At::pegtl_string,
                                store<ctr::param, ctr::beta, ctr::atwood>> >
               > {};

  //! geometry definition types
  struct geometry :
         sor< analytic_geometry,
              discrete_geometry > {};

  //! common to all physics
  struct physics_common :
         sor< process<kw::nstep::pegtl_string, store<ctr::incpar, ctr::nstep>>,
              process<kw::term::pegtl_string, store<ctr::incpar, ctr::term>>,
              process<kw::dt::pegtl_string, store<ctr::incpar, ctr::dt>>,
              process<kw::npar::pegtl_string, store<ctr::component, ctr::npar>>,
              process_quoted<kw::output::pegtl_string, set<ctr::io, ctr::physoutput>>,
              process<kw::pdfname::pegtl_string, set<ctr::io, ctr::pdf>>,
              process<kw::globname::pegtl_string, set<ctr::io, ctr::glob>>,
              process<kw::statname::pegtl_string, set<ctr::io, ctr::stats>>,
              process<kw::glbi::pegtl_string, store<ctr::interval, ctr::glob>>,
              process<kw::pdfi::pegtl_string, store<ctr::interval, ctr::pdf>>,
              process<kw::stai::pegtl_string, store<ctr::interval, ctr::plot>>,
              process<kw::ttyi::pegtl_string, store<ctr::interval, ctr::tty>>,
              process<kw::dmpi::pegtl_string, store<ctr::interval, ctr::dump>>
            > {};

  //! mass models
  struct mass :
         sor< beta > {};

  //! hydro models
  struct hydro :
         sor< slm > {};

  //! material mix models
  struct mix :
         sor< dir, gendir > {};

  //! turbulence frequency models
  struct freq :
         sor< freq_gamma > {};

  //! physics 'hommix' block
  struct hommix :
         ifmust< parse<kw::hommix::pegtl_string,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       mix,
                       statistics> > {};

  //! physics 'homrt' block
  struct homrt :
         ifmust< parse<kw::homrt::pegtl_string,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       mass,
                       hydro,
                       freq,
                       statistics> > {};

  //! physics 'homhydro' block
  struct homhydro :
         ifmust< parse<kw::homhydro::pegtl_string,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       hydro,
                       freq,
                       statistics> > {};

  //! physics 'spinsflow' block
  struct spinsflow :
         ifmust< parse<kw::spinsflow::pegtl_string,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       hydro,
                       freq,
                       mix> > {};

  //! physics types
  struct physics :
         sor< hommix,
              homhydro,
              homrt,
              spinsflow > {};

  //! main keywords
  struct keywords :
         sor< title,
              physics > {};

  //! ignore: comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore, unknown<Error::KEYWORD>> > {};

} // grm::
} // quinoa::

#endif // QuinoaGrammar_h
