//******************************************************************************
/*!
  \file      src/Control/QuinoaGrammar.h
  \author    J. Bakosi
  \date      Mon Sep  9 09:55:42 2013
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

#include <algorithm>

#include <Macro.h>
#include <QuinoaControlTypes.h>
#include <Option.h>
#include <Box.h>
#include <QuinoaKeywords.h>

namespace quinoa {

namespace grammar {

  using namespace pegtl;

  // State

  //! Bundle is where everything is stored during parsing
  using Stack = QuinoaControl;
  //! Out-of-struct storage of field ID for push_term
  static int field = 0;
  //! Geometry definition options
  static control::Option<select::Geometry> Geometry;
  //! Physics options
  static control::Option<select::Physics> Physics;
  //! Position options
  static control::Option<select::Position> Position;
  //! Mass options
  static control::Option<select::Mass> Mass;
  //! Hydro options
  static control::Option<select::Hydro> Hydro;
  //! Energy options
  static control::Option<select::Energy> Energy;
  //! Material mix options
  static control::Option<select::Mix> Mix;
  //! Turbulence frequency options
  static control::Option<select::Frequency> Frequency;
  //! Material mix rate options
  static control::Option<select::MixRate> MixRate;

  // Actions

  // convert string to 'type' via std::stringstream
  template< typename type >
  static type convert(const std::string& str) {
    std::stringstream ss(str);
    type num;
    ss >> num;
    return num;
  }

  // convert 'type' to string via std::stringstream
  template< typename type >
  static std::string convert(const type& val) {
    std::stringstream ss;
    ss << val;
    return ss.str();
  }

  // convert & store value in state at position given by tags
  template< typename... tags >
  struct cstore : action_base< cstore<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      // Convert to correct type and store element
      //stack.set<tags...>(
      //  convert< typename std::tuple_element<at,Stack>::type >(value);
    }
  };

  // store value in state at position given by tags (no conversion)
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.set<tags...>(value);
    }
  };

  // push value in state at position given by tags, converting 'to' type
  template< typename to, typename...tags >
  struct push : action_base< push<to,tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      // Convert to correct type and push element to vector
      //std::get<at>(stack).push_back(convert<to>(value));
    }
  };

  // start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value, Stack& stack) {
      //stack.get<control::stats>().push_back(std::vector<control::Term>());
      IGNORE(value);        // suppress compiler warning on unused variable
    }
  };

  // add matched value as Term into vector of Product in vector of statistics
  template< control::Quantity quantity, control::Moment moment, char name='\0' >
  struct push_term : action_base< push_term<quantity, moment, name> > {
    static void apply(const std::string& value, Stack& stack) {
// 
//       // if name is given, push name, otherwise push first char of value
//       char na(name ? name : value[0]);
//       // if name is given, it is triggered not user-requested
//       bool plot(name ? false : true);
//       // Use stats for shorthand of reference in bundle
//       std::vector<control::Product>& stats =
//         std::get<control::STATISTICS>(stack);
//       // Push term into current product
//       stats.back().push_back(control::Term(field, quantity, moment, na, plot));
// 
//       // If central moment, trigger mean
//       if (moment == control::Moment::CENTRAL) {
//         control::Term term(field,
//                            quantity,
//                            control::Moment::ORDINARY,
//                            toupper(na),
//                            false);
//         stats.insert(stats.end()-1, control::Product(1,term));
//       }
// 
//       field = 0;            // reset default field
    }
  };

  // save field ID so push_term can pick it up
  struct save_field : action_base< save_field > {
    static void apply(const std::string& value, Stack& stack) {
      field = convert<int>(value) - 1;  // numbering of field IDs start from 0
      IGNORE(stack);        // suppress compiler warning on unused variable
    }
  };

  template< class OptionType, typename...tags >
  struct store_option : action_base< store_option<OptionType, tags...> > {
    static void apply(const std::string& value, Stack& stack) {
//       control::Option<OptionType> Model;
//       // Issue warning if overwrite
//       if (boolstack[at]) {
//         std::cout << ">>> PARSER WARNING: Conflicting options defined in input "
//                      "file, overwriting \"" << Model.name(std::get<at>(stack))
//                   << "\" with \"" << Model.name(Model.value(value)) << "\""
//                   << std::endl;
//       }
//       std::get<at>(stack) = Model.value(value);
//       boolstack[at] = true;
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
  struct position : sor< keyword::pos_inviscid,
                         keyword::pos_viscous > {};

  // match all accepted as hydro models
  struct hydro : sor< keyword::hydro_slm,
                      keyword::hydro_glm > {};

  // match all accepted as mix models
  struct mix : sor< keyword::mix_iem,
                    keyword::mix_iecm,
                    keyword::mix_dir,
                    keyword::mix_gendir > {};

  // parse input padded by blank at left and space at right and if it matches
  // 'keywords', apply 'actions'
  template< class keywords, typename ... actions >
  struct parse :
         pad< ifapply< trim<keywords, space>, actions ... >, blank, space > {};

  // comment: start with '#' until eol
  struct comment :
         pad< trim<one<'#'>,eol>, blank, eol> {};

  // plow through block of 'tokens' until 'end' keyword
  template< typename ... tokens >
  struct block :
         until< read<keyword::end>, sor<comment, tokens ...> > {};

  // rng: one of the random number generators
  struct rng :
         sor< keyword::mkl_mcg31,
              keyword::mkl_r250,
              keyword::mkl_mrg32k3a,
              keyword::mkl_mcg59,
              keyword::mkl_wh,
              keyword::mkl_mt19937,
              keyword::mkl_mt2203,
              keyword::mkl_sfmt19937,
              keyword::mkl_sobol,
              keyword::mkl_niederr,
              keyword::mkl_iabstract,
              keyword::mkl_dabstract,
              keyword::mkl_sabstract,
              keyword::mkl_nondeterm > {};

  // number: optional sign followed by digits
  struct number :
         seq< opt< sor<one<'+'>, one<'-'>> >, digit> {};

  // plow through list of values between keywords 'key' and "end", calling
  // 'insert' for each if matches and allowing comments between values
  template< class key, class insert, class value = number >
  struct list :
         ifmust< read<key>,
                 until< read<keyword::end>, sor<comment, parse<value,insert>> >
               > {};

  // process 'keyword' and call its 'insert' action if matches 'keywords'
  template< class keyword, class insert, class keywords = alnum >
  struct process :
         ifmust< read<keyword>, parse<keywords,insert> > {};

  // moment: 'keyword' optionally followed by a digit, pushed to vector of terms
  template< class keyword, control::Quantity q, control::Moment m >
  struct moment :
         sor < ifapply<seq<keyword, ifapply<digit,save_field>>, push_term<q,m>>,
               ifapply<keyword, push_term<q,m>>
             > {};

  // terms recognized within an expectation and their mapping
  struct terms :
         sor< moment<keyword::transported_scalar,
                     control::Quantity::SCALAR,
                     control::Moment::ORDINARY>,
              moment<keyword::transported_scalar_fluctuation,
                     control::Quantity::SCALAR,
                     control::Moment::CENTRAL>,
              moment<keyword::velocity_x,
                     control::Quantity::VELOCITY_X,
                     control::Moment::ORDINARY>
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
                          store<control::title> >,
                 one<rbound>
               > {};

  // title
  struct title :
         ifmust< read<keyword::title>, parse_title<'"','"'> > {};

  // analytic_geometry block
  struct analytic_geometry:
         ifmust< parse< keyword::analytic_geometry,
                        store_option<select::Geometry,
                                     control::selected,
                                     control::geometry> > > {};

  // discrete_geometry block
  struct discrete_geometry:
         ifmust< parse< keyword::discrete_geometry,
                        store_option<select::Geometry,
                                     control::selected,
                                     control::geometry> > > {};

  // dir block
  struct dir :
         ifmust< parse< keyword::mix_dir,
                        store_option<select::Mix,
                                     control::selected,
                                     control::mix>>,
                 block< process<keyword::nscalar,
                                cstore<control::component, control::nscalar>>,
                        list<keyword::dir_B, push<control::param,
                                                  control::dirichlet,
                                                  control::b>>,
                        list<keyword::dir_S, push<control::param,
                                                  control::dirichlet,
                                                  control::S>>,
                        list<keyword::dir_kappa, push<control::param,
                                                      control::dirichlet,
                                                      control::kappa>> > > {};

  // gendir block
  struct gendir :
         ifmust< parse< keyword::mix_gendir,
                        store_option<select::Mix,
                                     control::selected,
                                     control::mix>>,
                 block< process<keyword::nscalar,
                                cstore<control::component, control::nscalar>>,
                        list<keyword::dir_B, push<control::param,
                                                  control::gendirichlet,
                                                  control::b>>,
                        list<keyword::dir_S, push<control::param,
                                                  control::gendirichlet,
                                                  control::S>>,
                        list<keyword::dir_kappa, push<control::param,
                                                      control::gendirichlet,
                                                      control::kappa>>,
                        list<keyword::gendir_C, push<control::param,
                                                     control::gendirichlet,
                                                     control::c>> > > {};

  // statistics block
  struct statistics :
         ifmust< read< keyword::statistics >,
                 block< parse_expectations<'<','>'> > > {};

  // slm block
  struct slm :
         ifmust< parse< keyword::hydro_slm,
                        store_option<select::Hydro,
                                     control::selected,
                                     control::hydro>,
                        // trigger estimating the diagonal of Reynolds-stress
                        start_product,
                        push_term<control::Quantity::VELOCITY_X,
                                  control::Moment::CENTRAL, 'u'>,
                        push_term<control::Quantity::VELOCITY_X,
                                  control::Moment::CENTRAL, 'u'>,
                        start_product,
                        push_term<control::Quantity::VELOCITY_Y,
                                  control::Moment::CENTRAL, 'v'>,
                        push_term<control::Quantity::VELOCITY_Y,
                                  control::Moment::CENTRAL, 'v'>,
                        start_product,
                        push_term<control::Quantity::VELOCITY_Z,
                                  control::Moment::CENTRAL, 'w'>,
                        push_term<control::Quantity::VELOCITY_Z,
                                  control::Moment::CENTRAL, 'w'>>,
                 block< process<keyword::SLM_C0, cstore<control::param,
                                                        control::slm,
                                                        control::c0>>,
                        process<keyword::nvelocity, cstore<control::component,
                                                           control::nvelocity>>
                      > > {};

  // freq_gamma block
  struct freq_gamma :
         ifmust< parse< keyword::freq_gamma,
                        store_option<select::Frequency,
                                     control::selected,
                                     control::frequency> >,
                 block< process<keyword::nfreq, cstore<control::component,
                                                       control::nfrequency>>,
                        process<keyword::freq_gamma_C1, cstore<control::param,
                                                               control::gamma,
                                                               control::c1>>,
                        process<keyword::freq_gamma_C2, cstore<control::param,
                                                               control::gamma,
                                                               control::c2>>,
                        process<keyword::freq_gamma_C3, cstore<control::param,
                                                               control::gamma,
                                                               control::c3>>,
                        process<keyword::freq_gamma_C4, cstore<control::param,
                                                               control::gamma,
                                                               control::c4>> >
               > {};

  // beta block
  struct beta :
         ifmust< parse< keyword::mass_beta,
                        store_option<select::Mass,
                                     control::selected,
                                     control::mass> >,
                 block< process<keyword::ndensity, cstore<control::component,
                                                          control::ndensity>>,
                        process<keyword::Beta_At, cstore<control::param,
                                                         control::beta,
                                                         control::atwood>> >
               > {};

  // geometry definition types
  struct geometry :
         sor< analytic_geometry,
              discrete_geometry > {};

  // common to all physics
  struct physics_common :
         sor< process<keyword::nstep, cstore<control::incpar, control::nstep>>,
              process<keyword::term, cstore<control::incpar, control::term>>,
              process<keyword::dt, cstore<control::incpar, control::dt>>,
              process<keyword::npar, cstore<control::component, control::npar>>,
              process<keyword::input, store<control::io, control::input>>,
              process<keyword::output, store<control::io, control::output>>,
              process<keyword::pdfname, store<control::io, control::pdf>>,
              process<keyword::globname, store<control::io, control::glob>>,
              process<keyword::statname, store<control::io, control::stats>>,
              process<keyword::glbi, cstore<control::interval, control::glob>>,
              process<keyword::pdfi, cstore<control::interval, control::pdf>>,
              process<keyword::stai, cstore<control::interval, control::stats>>,
              process<keyword::ttyi, cstore<control::interval, control::tty>>,
              process<keyword::dmpi, cstore<control::interval, control::dump>>
            > {};

  // hommix block
  struct hommix :
         ifmust< parse< keyword::hommix,
                        store_option<select::Physics,
                                     control::selected,
                                     control::physics> >,
                 block< geometry,
                        physics_common,
                        dir,
                        gendir,
                        statistics > > {};

  // homrt block
  struct homrt :
         ifmust< parse< keyword::homrt,
                        store_option<select::Physics,
                                     control::selected,
                                     control::physics> >,
                 block< geometry,
                        physics_common,
                        dir,
                        gendir,
                        slm,
                        beta,
                        statistics > > {};

  // homhydro block
  struct homhydro :
         ifmust< parse< keyword::homhydro,
                        store_option<select::Physics,
                                     control::selected,
                                     control::physics> >,
                 block< geometry,
                        physics_common,
                        slm,
                        freq_gamma,
                        statistics > > {};

  // spinsflow block
  struct spinsflow :
         ifmust< parse< keyword::spinsflow,
                        store_option<select::Physics,
                                     control::selected,
                                     control::physics> >,
                 block< geometry,
                        physics_common,
                        slm,
                        freq_gamma,
                        dir,
                        gendir,
                        beta > > {};

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

} // namespace grammar

} // namespace quinoa

#endif // QuinoaGrammar_h
