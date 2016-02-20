//******************************************************************************
/*!
  \file      src/Control/Breeze/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 02:20:28 PM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Breeze's input deck grammar definition
  \details   Breeze's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef BreezeInputDeckGrammar_h
#define BreezeInputDeckGrammar_h

#include "CommonGrammar.h"
#include "PEGTLParsed.h"
#include "Keywords.h"

#ifdef HAS_MKL
#include "MKLGrammar.h"
#endif

#include "RNGSSEGrammar.h"

namespace breeze {

extern ctr::InputDeck g_inputdeck_defaults;

//! Breeze input deck facilitating user input for computing fluid dynamics
namespace deck {

  //! \brief PEGTLParsed type specialized to Breeze's input deck parser
  //! \details PEGTLInputDeck is practically InputDeck equipped with PEGTL
  //!   location information so the location can be tracked during parsing.
  //! \author J. Bakosi
  using PEGTLInputDeck =
    tk::ctr::PEGTLParsed< ctr::InputDeck,
                          pegtl::file_input< ctr::Location >,
                          tag::cmd,
                          ctr::CmdLine >;

  //! \brief Specialization of tk::grm::use for Breeze's input deck parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3,
                            ctr::InputDeck::keywords4 >;

  // Breeze's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  // Breeze's InputDeck actions

  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults.
  //! \author J. Bakosi
  template< class Option, typename... tags >
  struct store_option : pegtl::action_base< store_option< Option, tags... > > {
    static void apply(const std::string& value, Stack& stack) {
      tk::grm::store_option< Stack, use, Option, ctr::InputDeck, tags... >
                           ( stack, value, g_inputdeck_defaults );
    }
  };

  // Breeze's InputDeck grammar

  //! \brief Scan selected option
  //! \author J. Bakosi
  template< typename keyword, typename option, typename... tags >
  struct select_option :
         tk::grm::scan< Stack,
                        typename keyword::pegtl_string,
                        store_option< option, tag::selected, tags... > >
         {};

  //! \brief Scan selected option and trigger action(s)
  //! \author J. Bakosi
  template< typename keyword, typename option, typename Tag,
            typename... triggers >
  struct select_option_and_trigger :
         tk::grm::scan< Stack, typename keyword::pegtl_string,
                        store_option< option, tag::selected, Tag >,
                        triggers... > {};

  //! \brief Scan and store MonteCarlo keyword and option
  //! \author J. Bakosi
  template< typename keyword >
  struct scan_montecarlo :
         select_option< keyword, ctr::MonteCarlo, tag::montecarlo > {};

  //! \brief Scan and store mass keyword and option
  //! \author J. Bakosi
  template< typename keyword >
  struct scan_mass :
         select_option< keyword, ctr::Mass, tag::mass > {};

  //! \brief Scan and store hydro keyword and option
  //! \author J. Bakosi
  template< typename keyword >
  struct scan_hydro :
         select_option< keyword, ctr::Hydro, tag::hydro > {};

  //! \brief Scan and store mix keyword and option
  //! \author J. Bakosi
  template< typename keyword >
  struct scan_mix :
         select_option< keyword, ctr::Mix, tag::mix > {};

  //! \brief Scan and store frequency keyword and option
  //! \author J. Bakosi
  template< typename keyword >
  struct scan_frequency :
         select_option< keyword, ctr::Frequency, tag::frequency > {};

  //! \brief Match and set fluctuating velocity in x direction
  //! \author J. Bakosi
  struct u :
         tk::grm::push_term< Stack, tk::ctr::Moment::CENTRAL, 'u' > {};

  //! \brief Match and set fluctuating velocity in y direction
  //! \author J. Bakosi
  struct v :
         tk::grm::push_term< Stack, tk::ctr::Moment::CENTRAL, 'v' > {};

  //! \brief Match and set fluctuating velocity in z direction
  //! \author J. Bakosi
  struct w :
         tk::grm::push_term< Stack, tk::ctr::Moment::CENTRAL, 'w' > {};

  //! \brief Match slm ... end block
  //! \author J. Bakosi
  struct slm :
         pegtl::ifmust<
           select_option_and_trigger<
             use< kw::hydro_slm >,
             ctr::Hydro,
             tag::hydro,
             // trigger Reynolds-stress diagonal
             tk::grm::start_vector< Stack, tag::stat >, u, u,
             tk::grm::start_vector< Stack, tag::stat >, v, v,
             tk::grm::start_vector< Stack, tag::stat >, w, w >,
           tk::grm::block<
              Stack,
              kw::end,
              tk::grm::parameter< Stack,
                                  kw::SLM_C0,
                                  pegtl::digit,
                                  tag::slm,
                                  tag::c0 >,
              tk::grm::component< Stack, kw::nvelocity, tag::hydro > > > {};

  //! Match common settings to all monte-carlo
  //! \author J. Bakosi
  struct montecarlo_common :
         pegtl::sor< tk::grm::discr< Stack, kw::npar, tag::npar >,
                     tk::grm::discr< Stack, kw::nstep, tag::nstep >,
                     tk::grm::discr< Stack, kw::term, tag::term >,
                     tk::grm::discr< Stack, kw::dt, tag::dt >,
                     tk::grm::interval< Stack, kw::glbi, tag::glob >,
                     tk::grm::interval< Stack, kw::ttyi, tag::tty >,
                     tk::grm::interval< Stack, kw::dmpi, tag::dump > > {};

  //! \brief Match the inside of rngs ... end block
  //! \author J. Bakosi
  struct rngs :
         pegtl::sor<
                     #ifdef HAS_MKL
                     tk::mkl::rngs< Stack, use,
                                    tag::selected, tag::rng,
                                    tag::param, tag::rngmkl >,
                     #endif
                     tk::rngsse::rngs< Stack, use,
                                       tag::selected, tag::rng,
                                       tag::param, tag::rngsse > > {};

  //! \brief Match mass models
  //! \author J. Bakosi
//   struct mass :
//          pegtl::sor< mass_beta > {};

  //! \brief Match hydro models
  //! \author J. Bakosi
  struct hydro :
         pegtl::sor< slm > {};

  //! \brief Match material mix models
  //! \author J. Bakosi
//   struct mix :
//          pegtl::sor< mix_dir, mix_gendir > {};

  //! \brief Match turbulence frequency models
  //! \author J. Bakosi
//   struct freq :
//          pegtl::sor< freq_gamma > {};

  //! \brief Match hommix ... end block
  //! \author J. Bakosi
  struct hommix :
         pegtl::ifmust< scan_montecarlo< kw::hommix >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        montecarlo_common,
                                        //mix,
                                        tk::grm::rngblock< Stack, use, rngs >,
                                        tk::grm::statistics< Stack, use > > > {};

  //! \brief Match homrt ... end block
  //! \author J. Bakosi
  struct homrt :
         pegtl::ifmust< scan_montecarlo< kw::homrt >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        montecarlo_common,
                                        //mass,
                                        hydro,
                                        //freq,
                                        tk::grm::rngblock< Stack, use, rngs >,
                                        tk::grm::statistics< Stack, use > > > {};

  //! \brief Match homhydro ... end block
  //! \author J. Bakosi
  struct homhydro :
         pegtl::ifmust< scan_montecarlo< kw::homhydro >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        montecarlo_common,
                                        hydro,
                                        //freq,
                                        tk::grm::rngblock< Stack, use, rngs >,
                                        tk::grm::statistics< Stack, use > > > {};

  //! \brief Match spinsflow ... end block
  //! \author J. Bakosi
  struct spinsflow :
         pegtl::ifmust< scan_montecarlo< kw::spinsflow >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        montecarlo_common,
                                        hydro,
                                        //freq,
                                        //mix,
                                        tk::grm::rngblock< Stack, use, rngs >,
                                        tk::grm::statistics< Stack, use > > > {};

  //! \brief Match all physics types
  //! \author J. Bakosi
  struct physics :
         pegtl::sor< hommix,
                     homhydro,
                     homrt,
                     spinsflow > {};

  //! \brief All keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use >, physics > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  //! \author J. Bakosi
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore > {};

} // deck::
} // breeze::

#endif // BreezeInputDeckGrammar_h
