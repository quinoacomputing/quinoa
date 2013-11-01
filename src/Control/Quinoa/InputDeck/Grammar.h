//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Thu Oct 31 07:23:48 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck grammar definition
  \details   Quinoa's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef QuinoaInputDeckGrammar_h
#define QuinoaInputDeckGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <Option.h>
#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Quinoa/Types.h>
#include <Quinoa/InputDeck/Keywords.h>

namespace quinoa {
namespace deck {

  using namespace pegtl;
  using namespace tk::grm;

  //! PEGTLParsed type specialized to Quinoa's input deck parser
  using PEGTLInputDeck = ctr::PEGTLParsed< ctr::InputDeck,
                                           file_input< ctr::Location >,
                                           ctr::cmd,
                                           ctr::CmdLine >;

  // Quinoa's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;
  //! Out-of-struct storage of field ID for pushing terms for statistics
  static int field = 0;

  // Quinoa's InputDeck actions

  //! start new product in vector of statistics
  struct start_product : action_base< start_product > {
    static void apply(const std::string& value, Stack& stack) {
      stack.push_back<ctr::stat>(ctr::Product());
      IGNORE(value);   // suppress compiler warning: parameter never referenced
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
      std::vector<ctr::Product>& stats = stack.get<ctr::stat>();
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
      tk::Option<OptionType> opt;
      //! Emit warning on overwrite
      if (stack.get<tags...>() != ctr::InputDeckDefaults.get<tags...>()) {
        std::cout << "\n>>> PARSER WARNING: Multiple definitions for '"
                  << opt.group() << "' option. Overwriting '"
                  << opt.name(stack.get<tags...>()) << "' with '"
                  << opt.name(opt.value(value)) << "'.\n\n";
      }
      stack.set<tags...>(opt.value(value));
    }
  };

  // Quinoa's InputDeck grammar

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
              unknown<Stack,Error::MOMENT>
            > {};

  //! plow through terms in expectation until character 'rbound'
  template< char rbound >
  struct expectation :
         until< one<rbound>, terms > {};

  //! plow through expectations between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_expectations :
         readkw< ifmust< one<lbound>, apply<start_product>, expectation<rbound> >
               > {};

  //! title
  struct title :
         ifmust< readkw<kw::title::pegtl_string>,
                 quoted<Stack,Set<Stack,ctr::title>> > {};

  //! analytic_geometry block
  struct analytic_geometry:
         ifmust< scan< kw::analytic_geometry::pegtl_string,
                       store_option<ctr::Geometry,
                                    ctr::selected,
                                    ctr::geometry> >,
                 block<Stack,kw::end::pegtl_string> > {};

  //! discrete_geometry block
  struct discrete_geometry:
         ifmust< scan< kw::discrete_geometry::pegtl_string,
                       store_option<ctr::Geometry,
                                    ctr::selected,
                                    ctr::geometry> >,
                 block<Stack,kw::end::pegtl_string> > {};

  //! dir block
  struct dir :
         ifmust< scan< kw::mix_dir::pegtl_string,
                       store_option<ctr::Mix, ctr::selected, ctr::mix> >,
                 block< Stack,
                        kw::end::pegtl_string,
                        process<Stack,
                                kw::nscalar::pegtl_string,
                                Store<Stack, ctr::component, ctr::nscalar>>,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_B::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::dirichlet,
                                           ctr::b> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_S::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::dirichlet,
                                           ctr::S> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_kappa::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::dirichlet,
                                           ctr::kappa> > > > {};

  //! gendir block
  struct gendir :
         ifmust< scan< kw::mix_gendir::pegtl_string,
                       store_option<ctr::Mix, ctr::selected, ctr::mix> >,
                 block< Stack,
                        kw::end::pegtl_string,
                        process<Stack,
                                kw::nscalar::pegtl_string,
                                Store<Stack, ctr::component, ctr::nscalar> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_B::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::gendirichlet,
                                           ctr::b> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_S::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::gendirichlet,
                                           ctr::S> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::dir_kappa::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::gendirichlet,
                                           ctr::kappa> >,
                        vector< Stack,
                                kw::end::pegtl_string,
                                kw::gendir_C::pegtl_string,
                                Store_back<Stack,
                                           ctr::param,
                                           ctr::gendirichlet,
                                           ctr::c>> > > {};

  //! statistics block
  struct statistics :
         ifmust< readkw<kw::statistics::pegtl_string>,
                 block<Stack,
                       kw::end::pegtl_string,
                       parse_expectations<'<','>'>> > {};

  //! slm block
  struct slm :
         ifmust< scan< kw::hydro_slm::pegtl_string,
                       store_option<ctr::Hydro, ctr::selected, ctr::hydro>,
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
                 block< Stack,
                        kw::end::pegtl_string,
                        process<Stack,
                                kw::SLM_C0::pegtl_string,
                                Store<Stack, ctr::param, ctr::slm, ctr::c0>>,
                        process<Stack,
                                kw::nvelocity::pegtl_string,
                                Store<Stack, ctr::component, ctr::nvelocity>>
                      > > {};

  //! freq_gamma block
  struct freq_gamma :
         ifmust< scan< kw::freq_gamma::pegtl_string,
                       store_option<ctr::Frequency,
                                    ctr::selected,
                                    ctr::frequency> >,
                 block< Stack,
                        kw::end::pegtl_string,
                        process<Stack,
                                kw::nfreq::pegtl_string,
                                Store<Stack, ctr::component, ctr::nfrequency>>,
                        process<Stack,
                                kw::freq_gamma_C1::pegtl_string,
                                Store<Stack, ctr::param, ctr::gamma, ctr::c1>>,
                        process<Stack,
                                kw::freq_gamma_C2::pegtl_string,
                                Store<Stack, ctr::param, ctr::gamma, ctr::c2>>,
                        process<Stack,
                                kw::freq_gamma_C3::pegtl_string,
                                Store<Stack, ctr::param, ctr::gamma, ctr::c3>>,
                        process<Stack,
                                kw::freq_gamma_C4::pegtl_string,
                                Store<Stack, ctr::param, ctr::gamma, ctr::c4>> >
               > {};

  //! beta block
  struct beta :
         ifmust< scan<kw::mass_beta::pegtl_string,
                      store_option<ctr::Mass, ctr::selected, ctr::mass>>,
                 block< Stack,
                        kw::end::pegtl_string,
                        process<Stack,
                                kw::ndensity::pegtl_string,
                                Store<Stack, ctr::component, ctr::ndensity>>,
                        process<Stack,
                                kw::Beta_At::pegtl_string,
                                Store<Stack,
                                      ctr::param,
                                      ctr::beta,
                                      ctr::atwood>> >
               > {};

  //! geometry definition types
  struct geometry :
         sor< analytic_geometry,
              discrete_geometry > {};

  //! common to all physics
  struct physics_common :
         sor< process<Stack,
                      kw::nstep::pegtl_string,
                      Store<Stack, ctr::incpar, ctr::nstep>>,
              process<Stack,
                      kw::term::pegtl_string,
                      Store<Stack, ctr::incpar, ctr::term>>,
              process<Stack,
                      kw::dt::pegtl_string,
                      Store<Stack, ctr::incpar, ctr::dt>>,
              process<Stack,
                      kw::npar::pegtl_string,
                      Store<Stack, ctr::component, ctr::npar>>,
              process<Stack,
                      kw::glbi::pegtl_string,
                      Store<Stack, ctr::interval, ctr::glob>>,
              process<Stack,
                      kw::pdfi::pegtl_string,
                      Store<Stack, ctr::interval, ctr::pdf>>,
              process<Stack,
                      kw::stai::pegtl_string,
                      Store<Stack, ctr::interval, ctr::plot>>,
              process<Stack,
                      kw::ttyi::pegtl_string,
                      Store<Stack, ctr::interval, ctr::tty>>,
              process<Stack,
                      kw::dmpi::pegtl_string,
                      Store<Stack, ctr::interval, ctr::dump>>
            > {};

  //! rng block
  struct rng :
         ifmust< scan< kw::mkl_mcg31::pegtl_string,
                       store_option<ctr::RNG, ctr::selected, ctr::rng>>,
                 block< Stack,
                        kw::end::pegtl_string,
                        process< Stack,
                                 kw::seed::pegtl_string,
                                 Store< Stack, ctr::param,
                                               ctr::rng,
                                               ctr::seed> > > > {};

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
         ifmust< scan<kw::hommix::pegtl_string,
                      store_option<ctr::Physics, ctr::selected, ctr::physics>>,
                 block<Stack,
                       kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       mix,
                       rng,
                       statistics> > {};

  //! physics 'homrt' block
  struct homrt :
         ifmust< scan<kw::homrt::pegtl_string,
                      store_option<ctr::Physics, ctr::selected, ctr::physics>>,
                 block<Stack,
                       kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       mass,
                       hydro,
                       freq,
                       statistics> > {};

  //! physics 'homhydro' block
  struct homhydro :
         ifmust< scan<kw::homhydro::pegtl_string,
                      store_option<ctr::Physics, ctr::selected, ctr::physics>>,
                 block<Stack,
                       kw::end::pegtl_string,
                       geometry,
                       physics_common,
                       hydro,
                       freq,
                       statistics> > {};

  //! physics 'spinsflow' block
  struct spinsflow :
         ifmust< scan<kw::spinsflow::pegtl_string,
                      store_option<ctr::Physics, ctr::selected, ctr::physics>>,
                 block<Stack,
                       kw::end::pegtl_string,
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
         until< eof, sor<keywords, ignore, unknown<Stack,Error::KEYWORD>> > {};

} // deck::
} // quinoa::

#endif // QuinoaInputDeckGrammar_h
