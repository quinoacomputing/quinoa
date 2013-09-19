//******************************************************************************
/*!
  \file      src/Control/QuinoaGrammar.h
  \author    J. Bakosi
  \date      Thu Sep 19 10:36:22 2013
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
#include <QuinoaControlTypes.h>
#include <Option.h>
#include <QuinoaKeywords.h>

namespace quinoa {
//! Grammar definition: state, actions, grammar definition
namespace grm {

  using namespace pegtl;

  // State

  //! Bundle is where everything is stored during parsing
  using Stack = QuinoaControl;
  //! Out-of-struct storage of field ID for pushing terms for statistics
  static int field = 0;

  // Actions

  // put value in state at position given by tags without conversion
  template< typename... tags >
  struct set : action_base< set<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.set<tags...>(value);
    }
  };

  // convert and put value in state at position given by tags
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store<tags...>(value);
    }
  };

  // convert and push back value to vector in state at position given by tags
  template< typename...tags >
  struct store_back : action_base< store_back<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store_back<tags...>(value);
    }
  };

  // start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value, Stack& stack) {
      stack.push_back<ctr::stats>(ctr::Product());
      IGNORE(value);        // suppress compiler warning on unused variable
    }
  };

  // add matched value as Term into vector of Product in vector of statistics
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

  // save field ID so push_term can pick it up
  struct save_field : action_base< save_field > {
    static void apply(const std::string& value, Stack& stack) {
      field = stack.convert<int>(value) - 1;  // field ID numbering start from 0
    }
  };

  // convert and put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : action_base< store_option<OptionType, tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      ctr::Option<OptionType> Model;
      stack.set<tags...>(Model.value(value));
    }
  };

  // Grammar

  // read 'token' until 'erased' trimming 'erased'
  template< class token, class erased >
  struct trim :
         seq< token, until< at<erased> > > {};

  // read 'token' padded by blank at left and space at right
  template< class token >
  struct read :
         pad< trim<token, space>, blank, space > {};

  // match all accepted as position models
  struct position : sor< kw::pos_inviscid,
                         kw::pos_viscous > {};

  // match all accepted as hydro models
  struct hydro : sor< kw::hydro_slm,
                      kw::hydro_glm > {};

  // match all accepted as mix models
  struct mix : sor< kw::mix_iem,
                    kw::mix_iecm,
                    kw::mix_dir,
                    kw::mix_gendir > {};

  // parse input padded by blank at left and space at right and if it matches
  // 'keywords', apply 'actions'
  template< class keywords, typename... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions... >, blank, space > {};

  // comment: start with '#' until eol
  struct comment :
         pad< trim<one<'#'>,eol>, blank, eol> {};

  // plow through block of 'tokens' until 'end' keyword
  template< typename... tokens >
  struct block :
         until< read<kw::end>, sor<comment, tokens...> > {};

  // rng: one of the random number generators
  struct rng :
         sor< kw::mkl_mcg31,
              kw::mkl_r250,
              kw::mkl_mrg32k3a,
              kw::mkl_mcg59,
              kw::mkl_wh,
              kw::mkl_mt19937,
              kw::mkl_mt2203,
              kw::mkl_sfmt19937,
              kw::mkl_sobol,
              kw::mkl_niederr,
              kw::mkl_iabstract,
              kw::mkl_dabstract,
              kw::mkl_sabstract,
              kw::mkl_nondeterm > {};

  // number: optional sign followed by digits
  struct number :
         seq< opt< sor<one<'+'>, one<'-'>> >, digit> {};

  // plow through list of values between keywords 'key' and "end", calling
  // 'insert' for each if matches and allowing comments between values
  template< class key, class insert, class value = number >
  struct list :
         ifmust< read<key>,
                 until< read<kw::end>, sor<comment, parse<value,insert>> >
               > {};

  // process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<keywords,insert> > {};

  // moment: 'keyword' optionally followed by a digit, pushed to vector of terms
  template< class keyword, ctr::Quantity q, ctr::Moment m >
  struct moment :
         sor < ifapply<seq<keyword, ifapply<digit,save_field>>, push_term<q,m>>,
               ifapply<keyword, push_term<q,m>> > {};

  // terms recognized within an expectation and their mapping
  struct terms :
         sor< moment<kw::transported_scalar,
                     ctr::Quantity::SCALAR,
                     ctr::Moment::ORDINARY>,
              moment<kw::transported_scalar_fluctuation,
                     ctr::Quantity::SCALAR,
                     ctr::Moment::CENTRAL>,
              moment<kw::velocity_x,
                     ctr::Quantity::VELOCITY_X,
                     ctr::Moment::ORDINARY>
            > {};

  // plow through terms in expectation until character 'rbound'
  template< char rbound >
  struct expectation :
         until< one<rbound>, terms > {};

  // plow through expectations between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_expectations :
         read< ifmust< one<lbound>, apply<start_product>, expectation<rbound> >
             > {};

  // scan problem title between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_title :
         ifmust< one<lbound>,
                 ifapply< trim<not_one<lbound>, one<rbound>>,
                          set<ctr::title> >,
                 one<rbound>
               > {};

  // title
  struct title :
         ifmust< read<kw::title>, parse_title<'"','"'> > {};

  // analytic_geometry block
  struct analytic_geometry:
         ifmust< parse< kw::analytic_geometry,
                        store_option<sel::Geometry,
                                     ctr::selected,
                                     ctr::geometry> >, block<> > {};

  // discrete_geometry block
  struct discrete_geometry:
         ifmust< parse< kw::discrete_geometry,
                        store_option<sel::Geometry,
                                     ctr::selected,
                                     ctr::geometry> >, block<> > {};

  // dir block
  struct dir :
         ifmust< parse< kw::mix_dir,
                        store_option<sel::Mix, ctr::selected, ctr::mix> >,
                 block< process<kw::nscalar,
                                store<ctr::component, ctr::nscalar>>,
                        list< kw::dir_B,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::b> >,
                        list< kw::dir_S,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::S> >,
                        list< kw::dir_kappa,
                              store_back<ctr::param,
                                         ctr::dirichlet,
                                         ctr::kappa> > > > {};

  // gendir block
  struct gendir :
         ifmust< parse< kw::mix_gendir,
                        store_option<sel::Mix, ctr::selected, ctr::mix> >,
                 block< process< kw::nscalar,
                                 store<ctr::component, ctr::nscalar> >,
                        list< kw::dir_B,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::b> >,
                        list< kw::dir_S,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::S> >,
                        list< kw::dir_kappa,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::kappa> >,
                        list< kw::gendir_C,
                              store_back<ctr::param,
                                         ctr::gendirichlet,
                                         ctr::c>> > > {};

  // statistics block
  struct statistics :
         ifmust< read< kw::statistics >,
                 block< parse_expectations<'<','>'> > > {};

  // slm block
  struct slm :
         ifmust< parse< kw::hydro_slm,
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
                 block< process<kw::SLM_C0,
                                store<ctr::param, ctr::slm, ctr::c0>>,
                        process<kw::nvelocity,
                                store<ctr::component, ctr::nvelocity>>
                      > > {};

  // freq_gamma block
  struct freq_gamma :
         ifmust< parse< kw::freq_gamma, store_option<sel::Frequency,
                                                     ctr::selected,
                                                     ctr::frequency> >,
                 block< process<kw::nfreq,
                                store<ctr::component, ctr::nfrequency>>,
                        process<kw::freq_gamma_C1,
                                store<ctr::param, ctr::gamma, ctr::c1>>,
                        process<kw::freq_gamma_C2,
                                store<ctr::param, ctr::gamma, ctr::c2>>,
                        process<kw::freq_gamma_C3,
                                store<ctr::param, ctr::gamma, ctr::c3>>,
                        process<kw::freq_gamma_C4,
                                store<ctr::param, ctr::gamma, ctr::c4>> >
               > {};

  // beta block
  struct beta :
         ifmust< parse<kw::mass_beta,
                       store_option<sel::Mass, ctr::selected, ctr::mass>>,
                 block< process<kw::ndensity,
                                store<ctr::component, ctr::ndensity>>,
                        process<kw::Beta_At,
                                store<ctr::param, ctr::beta, ctr::atwood>> >
               > {};

  // geometry definition types
  struct geometry :
         sor< analytic_geometry,
              discrete_geometry > {};

  // common to all physics
  struct physics_common :
         sor< process<kw::nstep, store<ctr::incpar, ctr::nstep>>,
              process<kw::term, store<ctr::incpar, ctr::term>>,
              process<kw::dt, store<ctr::incpar, ctr::dt>>,
              process<kw::npar, store<ctr::component, ctr::npar>>,
              process<kw::input, set<ctr::io, ctr::input>>,
              process<kw::output, set<ctr::io, ctr::output>>,
              process<kw::pdfname, set<ctr::io, ctr::pdf>>,
              process<kw::globname, set<ctr::io, ctr::glob>>,
              process<kw::statname, set<ctr::io, ctr::stats>>,
              process<kw::glbi, store<ctr::interval, ctr::glob>>,
              process<kw::pdfi, store<ctr::interval, ctr::pdf>>,
              process<kw::stai, store<ctr::interval, ctr::plot>>,
              process<kw::ttyi, store<ctr::interval, ctr::tty>>,
              process<kw::dmpi, store<ctr::interval, ctr::dump>>
            > {};

  // hommix block
  struct hommix :
         ifmust< parse<kw::hommix,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block< geometry,
                        physics_common,
                        dir,
                        gendir,
                        statistics > > {};

  // homrt block
  struct homrt :
         ifmust< parse<kw::homrt,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<geometry,
                       physics_common,
                       dir,
                       gendir,
                       slm,
                       beta,
                       statistics> > {};

  // homhydro block
  struct homhydro :
         ifmust< parse<kw::homhydro,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<geometry,
                       physics_common,
                       slm,
                       freq_gamma,
                       statistics> > {};

  // spinsflow block
  struct spinsflow :
         ifmust< parse<kw::spinsflow,
                       store_option<sel::Physics, ctr::selected, ctr::physics>>,
                 block<geometry,
                       physics_common,
                       slm,
                       freq_gamma,
                       dir,
                       gendir,
                       beta> > {};

  // physics
  struct physics :
         sor< hommix,
              homhydro,
              homrt,
              spinsflow > {};

  // main keywords
  struct keywords :
         sor< title,
              physics > {};

  // ignore: comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  // entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore> > {};

} // grm::
} // quinoa::

#endif // QuinoaGrammar_h
