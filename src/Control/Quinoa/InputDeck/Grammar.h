//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Mon 11 Nov 2013 10:34:37 AM MST
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

  //! PEGTLParsed type specialized to Quinoa's input deck parser
  using PEGTLInputDeck = ctr::PEGTLParsed< ctr::InputDeck,
                                           pegtl::file_input< ctr::Location >,
                                           ctr::cmd,
                                           ctr::CmdLine >;

  // Quinoa's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;
  //! Out-of-struct storage of field ID for pushing terms for statistics
  static int field = 0;

  // Quinoa's InputDeck actions

  //! start new product in vector of statistics
  struct start_product : pegtl::action_base< start_product > {
    static void apply(const std::string& value, Stack& stack) {
      stack.push_back<ctr::stat>(ctr::Product());
      IGNORE(value);   // suppress compiler warning: parameter never referenced
    }
  };

  //! add matched value as Term into vector of Product in vector of statistics
  template< ctr::Quantity q, ctr::Moment m, char name='\0' >
  struct push_term : pegtl::action_base< push_term<q, m, name> > {
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
  struct save_field : pegtl::action_base< save_field > {
    static void apply(const std::string& value, Stack& stack) {
      field = stack.convert<int>(value) - 1;  // field ID numbering start from 0
    }
  };

  //! put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : pegtl::action_base< store_option<OptionType,tags...> > {
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
         pegtl::sor < pegtl::ifapply< pegtl::seq< keyword,
                                                  pegtl::ifapply< pegtl::digit,
                                                                  save_field > >,
                                      push_term< q, m > >,
                      pegtl::ifapply< keyword, push_term< q, m > > > {};

  //! terms recognized within an expectation and their mapping
  struct terms :
         pegtl::sor< moment< kw::transported_scalar::pegtl_string,
                             ctr::Quantity::SCALAR,
                             ctr::Moment::ORDINARY >,
                     moment< kw::transported_scalar_fluctuation::pegtl_string,
                             ctr::Quantity::SCALAR,
                             ctr::Moment::CENTRAL >,
                     moment< kw::velocity_x::pegtl_string,
                             ctr::Quantity::VELOCITY_X,
                             ctr::Moment::ORDINARY >,
                     tk::grm::unknown< Stack, tk::grm::Error::MOMENT > > {};

  //! plow through terms in expectation until character 'rbound'
  template< char rbound >
  struct expectation :
         pegtl::until< pegtl::one< rbound >, terms > {};

  //! plow through expectations between characters 'lbound' and 'rbound'
  template< char lbound, char rbound >
  struct parse_expectations :
         tk::grm::readkw< pegtl::ifmust< pegtl::one< lbound >,
                                         pegtl::apply< start_product >,
                                         expectation< rbound > > > {};

  //! title
  struct title :
         pegtl::ifmust< tk::grm::readkw< tk::kw::title::pegtl_string >,
                                         tk::grm::quoted<
                                           Stack,
                                           tk::grm::Set< Stack,
                                                         ctr::title > > > {};

  //! analytic_geometry block
  struct analytic_geometry:
         pegtl::ifmust< tk::grm::scan< kw::analytic_geometry::pegtl_string,
                                       store_option< ctr::Geometry,
                                                     ctr::selected,
                                                     ctr::geometry > >,
                        tk::grm::block< Stack > > {};

  //! discrete_geometry block
  struct discrete_geometry:
         pegtl::ifmust< tk::grm::scan< kw::discrete_geometry::pegtl_string,
                                       store_option< ctr::Geometry,
                                                     ctr::selected,
                                                     ctr::geometry > >,
                        tk::grm::block< Stack > > {};

  //! dir block
  struct dir :
         pegtl::ifmust< tk::grm::scan< kw::mix_dir::pegtl_string,
                                       store_option< ctr::Mix,
                                                     ctr::selected,
                                                     ctr::mix > >,
                        tk::grm::block< Stack,
                                        tk::grm::process<
                                          Stack,
                                          kw::nscalar::pegtl_string,
                                          tk::grm::Store< Stack,
                                                          ctr::component,
                                                          ctr::nscalar > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_B::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::dirichlet,
                                                               ctr::b > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_S::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::dirichlet,
                                                               ctr::S > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_kappa::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::dirichlet,
                                                               ctr::kappa > > > > {};

  //! gendir block
  struct gendir :
         pegtl::ifmust< tk::grm::scan< kw::mix_gendir::pegtl_string,
                                       store_option< ctr::Mix,
                                                     ctr::selected,
                                                     ctr::mix > >,
                        tk::grm::block< Stack,
                                        tk::grm::process<
                                          Stack,
                                          kw::nscalar::pegtl_string,
                                          tk::grm::Store< Stack,
                                                          ctr::component,
                                                          ctr::nscalar > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_B::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::gendirichlet,
                                                               ctr::b > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_S::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::gendirichlet,
                                                               ctr::S > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::dir_kappa::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::gendirichlet,
                                                               ctr::kappa > >,
                                        tk::grm::vector<
                                          Stack,
                                          kw::gendir_C::pegtl_string,
                                          tk::grm::Store_back< Stack,
                                                               ctr::param,
                                                               ctr::gendirichlet,
                                                               ctr::c > > > > {};

  //! statistics block
  struct statistics :
         pegtl::ifmust< tk::grm::readkw< kw::statistics::pegtl_string >,
                        tk::grm::block< Stack, parse_expectations<'<','>'> > > {};

  //! slm block
  struct slm :
         pegtl::ifmust< tk::grm::scan<
                          kw::hydro_slm::pegtl_string,
                          store_option< ctr::Hydro,
                                        ctr::selected,
                                        ctr::hydro >,
                          // trigger estimating the diagonal of Reynolds-stress
                          start_product,
                          push_term< ctr::Quantity::VELOCITY_X,
                                     ctr::Moment::CENTRAL, 'u' >,
                          push_term< ctr::Quantity::VELOCITY_X,
                                     ctr::Moment::CENTRAL, 'u' >,
                          start_product,
                          push_term< ctr::Quantity::VELOCITY_Y,
                                     ctr::Moment::CENTRAL, 'v' >,
                          push_term< ctr::Quantity::VELOCITY_Y,
                                     ctr::Moment::CENTRAL, 'v' >,
                          start_product,
                          push_term< ctr::Quantity::VELOCITY_Z,
                                     ctr::Moment::CENTRAL, 'w' >,
                          push_term< ctr::Quantity::VELOCITY_Z,
                                     ctr::Moment::CENTRAL, 'w'> >,
                          tk::grm::block<
                            Stack,
                            tk::grm::process< Stack,
                                              kw::SLM_C0::pegtl_string,
                                              tk::grm::Store< Stack,
                                                              ctr::param,
                                                              ctr::slm,
                                                              ctr::c0 > >,
                            tk::grm::process< Stack,
                                              kw::nvelocity::pegtl_string,
                                              tk::grm::Store< Stack,
                                                              ctr::component,
                                                              ctr::nvelocity > > > > {};

  //! freq_gamma block
  struct freq_gamma :
         pegtl::ifmust< tk::grm::scan< kw::freq_gamma::pegtl_string,
                                       store_option< ctr::Frequency,
                                                     ctr::selected,
                                                     ctr::frequency > >,
                        tk::grm::block<
                          Stack,
                          tk::grm::process<
                            Stack,
                            kw::nfreq::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::component,
                                            ctr::nfrequency > >,
                          tk::grm::process<
                            Stack,
                            kw::freq_gamma_C1::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::param,
                                            ctr::gamma,
                                            ctr::c1 > >,
                          tk::grm::process<
                            Stack,
                            kw::freq_gamma_C2::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::param,
                                            ctr::gamma,
                                            ctr::c2 > >,
                          tk::grm::process<
                            Stack,
                            kw::freq_gamma_C3::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::param,
                                            ctr::gamma,
                                            ctr::c3 > >,
                          tk::grm::process<
                            Stack,
                            kw::freq_gamma_C4::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::param,
                                            ctr::gamma,
                                            ctr::c4 > > > > {};

  //! beta block
  struct beta :
         pegtl::ifmust< tk::grm::scan< kw::mass_beta::pegtl_string,
                                       store_option< ctr::Mass,
                                                     ctr::selected,
                                                     ctr::mass > >,
                        tk::grm::block<
                          Stack,
                          tk::grm::process<
                            Stack,
                            kw::ndensity::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::component,
                                            ctr::ndensity > >,
                          tk::grm::process<
                            Stack,
                            kw::Beta_At::pegtl_string,
                            tk::grm::Store< Stack,
                                            ctr::param,
                                            ctr::beta,
                                            ctr::atwood > > > > {};

  //! geometry definition types
  struct geometry :
         pegtl::sor< analytic_geometry,
                     discrete_geometry > {};

  //! common to all physics
  struct physics_common :
         pegtl::sor< tk::grm::process< Stack,
                                       kw::nstep::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::incpar,
                                                       ctr::nstep > >,
                     tk::grm::process< Stack,
                                       kw::term::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::incpar,
                                                       ctr::term > >,
                     tk::grm::process< Stack,
                                       kw::dt::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::incpar,
                                                       ctr::dt > >,
                     tk::grm::process< Stack,
                                       kw::npar::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::component,
                                                       ctr::npar > >,
                     tk::grm::process< Stack,
                                       kw::glbi::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::interval,
                                                       ctr::glob > >,
                     tk::grm::process< Stack,
                                       kw::pdfi::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::interval,
                                                       ctr::pdf > >,
                     tk::grm::process< Stack,
                                       kw::stai::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::interval,
                                                       ctr::plot > >,
                     tk::grm::process< Stack,
                                       kw::ttyi::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::interval,
                                                       ctr::tty > >,
                     tk::grm::process< Stack,
                                       kw::dmpi::pegtl_string,
                                       tk::grm::Store< Stack,
                                                       ctr::interval,
                                                       ctr::dump > > > {};

  //! mklrngs block
  struct mklrngs :
         pegtl::ifmust< tk::grm::scan< tk::grm::mklrng,
                                       store_option< ctr::RNG,
                                                     ctr::selected,
                                                     ctr::rng > >,
                        tk::grm::block<
                          Stack,
                          tk::grm::process< Stack,
                                            tk::kw::seed::pegtl_string,
                                            tk::grm::Store< Stack,
                                                            ctr::param,
                                                            ctr::rng,
                                                            ctr::seed > > > > {};

  //! rngs
  struct rngs :
         pegtl::sor< mklrngs > {};

  //! mass models
  struct mass :
         pegtl::sor< beta > {};

  //! hydro models
  struct hydro :
         pegtl::sor< slm > {};

  //! material mix models
  struct mix :
         pegtl::sor< dir, gendir > {};

  //! turbulence frequency models
  struct freq :
         pegtl::sor< freq_gamma > {};

  //! scan and store physics keyword and option
  template< typename keyword >
  struct scan_physics :
         tk::grm::scan< keyword,
                        store_option< ctr::Physics,
                                      ctr::selected,
                                      ctr::physics > > {};

  //! physics 'hommix' block
  struct hommix :
         pegtl::ifmust< scan_physics< kw::hommix::pegtl_string >,
                        tk::grm::block< Stack,
                                        geometry,
                                        physics_common,
                                        mix,
                                        rngs,
                                        statistics > > {};

  //! physics 'homrt' block
  struct homrt :
         pegtl::ifmust< scan_physics< kw::homrt::pegtl_string >,
                        tk::grm::block< Stack,
                                        geometry,
                                        physics_common,
                                        mass,
                                        hydro,
                                        freq,
                                        rngs,
                                        statistics > > {};

  //! physics 'homhydro' block
  struct homhydro :
         pegtl::ifmust< scan_physics< kw::homhydro::pegtl_string >,
                        tk::grm::block< Stack,
                                        geometry,
                                        physics_common,
                                        hydro,
                                        freq,
                                        rngs,
                                        statistics > > {};

  //! physics 'spinsflow' block
  struct spinsflow :
         pegtl::ifmust< scan_physics< kw::spinsflow::pegtl_string >,
                        tk::grm::block< Stack,
                                        geometry,
                                        physics_common,
                                        hydro,
                                        freq,
                                        mix,
                                        rngs,
                                        statistics > > {};

  //! physics types
  struct physics :
         pegtl::sor< hommix,
                     homhydro,
                     homrt,
                     spinsflow > {};

  //! main keywords
  struct keywords :
         pegtl::sor< title, physics > {};

  //! ignore: comments and empty lines
  struct ignore :
         pegtl::sor< tk::grm::comment,
                     pegtl::until< pegtl::eol, pegtl::space > > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< Stack, keywords, ignore > {};

} // deck::
} // quinoa::

#endif // QuinoaInputDeckGrammar_h
