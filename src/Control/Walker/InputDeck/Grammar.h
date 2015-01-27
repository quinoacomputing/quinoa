//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Fri 23 Jan 2015 06:35:37 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's input deck grammar definition
  \details   Walker's input deck grammar definition. We use the [Parsing
    Expression Grammar Template Library (PEGTL)]
    (https://code.google.com/p/pegtl/wiki/PEGTL0) to create the grammar and the
    associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL.
    Word of advice: read from the bottom up.
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

  //! \brief Specialization of tk::grm::use for Walker's input deck parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3,
                            ctr::InputDeck::keywords4,
                            ctr::InputDeck::keywords5 >;

  // Walker's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;
 
  // Walker's InputDeck actions

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

  // Walker's InputDeck grammar

  //! scan and store_back sde keyword and option
  template< typename keyword >
  struct scan_sde :
         tk::grm::scan< Stack, typename keyword::pegtl_string,
                        tk::grm::store_back_option< Stack,
                                                    use,
                                                    ctr::DiffEq,
                                                    tag::selected,
                                                    tag::diffeq > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discr< Stack, use< kw::npar >, tag::npar >,
                     tk::grm::discr< Stack, use< kw::nstep >, tag::nstep >,
                     tk::grm::discr< Stack, use< kw::term >, tag::term >,
                     tk::grm::discr< Stack, use< kw::dt >, tag::dt >,
                     tk::grm::interval< Stack, use< kw::ttyi >, tag::tty > > {};

  //! rngs
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

  //! Diagonal Ornstein-Uhlenbeck SDE
  struct diag_ou :
         pegtl::ifmust< scan_sde< use< kw::diag_ou > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::diagou,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::diagou >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::diagou,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::diagou,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::diagou,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_sigma >,
                                          tag::diagou,
                                          tag::sigma >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_theta >,
                                          tag::diagou,
                                          tag::theta >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_mu >,
                                          tag::diagou,
                                          tag::mu > > > {};

  //! Ornstein-Uhlenbeck SDE
  struct ornstein_uhlenbeck :
         pegtl::ifmust< scan_sde< use< kw::ornstein_uhlenbeck > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::ou,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::ou >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::ou,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::ou,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::ou,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_sigma >,
                                          tag::ou,
                                          tag::sigma >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_theta >,
                                          tag::ou,
                                          tag::theta >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_mu >,
                                          tag::ou,
                                          tag::mu > > > {};

  //! Skew-normal SDE
  struct skewnormal :
         pegtl::ifmust< scan_sde< use< kw::skewnormal > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::skewnormal,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::skewnormal >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::skewnormal,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::skewnormal,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::skewnormal,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_T >,
                                          tag::skewnormal,
                                          tag::timescale >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_sigma >,
                                          tag::skewnormal,
                                          tag::sigma >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_lambda >,
                                          tag::skewnormal,
                                          tag::lambda > > > {};

  //! Beta SDE
  struct beta :
         pegtl::ifmust<
           scan_sde< use< kw::beta > >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::depvar< Stack,
                                            use,
                                            tag::beta,
                                            tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::beta >,
                           tk::grm::rng< Stack,
                                         use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::beta,
                                         tag::rng >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::init >,
                                            tk::ctr::InitPolicy,
                                            tag::beta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            tk::ctr::CoeffPolicy,
                                            tag::beta,
                                            tag::coeffpolicy >,
                           tk::grm::parameter_vector<
                             Stack,
                             use,
                             use< kw::sde_b >,
                             tag::beta,
                             tag::b >,
                           tk::grm::parameter_vector<
                             Stack,
                             use,
                             use< kw::sde_S >,
                             tag::beta,
                             tag::S >,
                           tk::grm::parameter_vector<
                             Stack,
                             use,
                             use< kw::sde_kappa >,
                             tag::beta,
                             tag::kappa > > > {};

  //! Gamma SDE
  struct gamma :
         pegtl::ifmust< scan_sde< use< kw::gamma > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::gamma,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::gamma >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::gamma,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::gamma,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::gamma,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_b >,
                                          tag::gamma,
                                          tag::b >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_S >,
                                          tag::gamma,
                                          tag::S >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_kappa >,
                                          tag::gamma,
                                          tag::kappa > > > {};

  //! Dirichlet SDE
  struct dirichlet :
         pegtl::ifmust< scan_sde< use< kw::dirichlet > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::dirichlet,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::dirichlet >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::dirichlet,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::dirichlet,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::dirichlet,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_b >,
                                          tag::dirichlet,
                                          tag::b >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_S >,
                                          tag::dirichlet,
                                          tag::S >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_kappa >,
                                          tag::dirichlet,
                                          tag::kappa > > > {};
  //! Generalized Dirichlet SDE
  struct gendir :
         pegtl::ifmust< scan_sde< use< kw::gendir > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::gendir,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::gendir >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::gendir,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::gendir,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::gendir,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_b >,
                                          tag::gendir,
                                          tag::b >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_S >,
                                          tag::gendir,
                                          tag::S >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_kappa >,
                                          tag::gendir,
                                          tag::kappa >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_c >,
                                          tag::gendir,
                                          tag::c > > > {};

  //! Wright-Fisher SDE
  struct wright_fisher :
         pegtl::ifmust< scan_sde< use< kw::wrightfisher > >,
                        tk::grm::block< Stack,
                                        use< kw::end >,
                                        tk::grm::depvar< Stack,
                                                         use,
                                                         tag::wrightfisher,
                                                         tag::depvar >,
                                        tk::grm::component< Stack,
                                                            use< kw::ncomp >,
                                                            tag::wrightfisher >,
                                        tk::grm::rng< Stack,
                                                      use,
                                                      use< kw::rng >,
                                                      tk::ctr::RNG,
                                                      tag::wrightfisher,
                                                      tag::rng >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::init >,
                                                         tk::ctr::InitPolicy,
                                                         tag::wrightfisher,
                                                         tag::initpolicy >,
                                        tk::grm::policy< Stack,
                                                         use,
                                                         use< kw::coeff >,
                                                         tk::ctr::CoeffPolicy,
                                                         tag::wrightfisher,
                                                         tag::coeffpolicy >,
                                        tk::grm::parameter_vector<
                                          Stack,
                                          use,
                                          use< kw::sde_omega >,
                                          tag::wrightfisher,
                                          tag::omega > > > {};

  //! stochastic differential equations
  struct sde :
         pegtl::sor< dirichlet,
                     gendir,
                     wright_fisher,
                     ornstein_uhlenbeck,
                     diag_ou,
                     skewnormal,
                     gamma,
                     beta > {};

  //! 'walker' block
  struct walker :
         pegtl::ifmust<
           tk::grm::readkw< use< kw::walker >::pegtl_string >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           discretization_parameters,
                           sde,
                           tk::grm::rngblock< Stack, use, rngs >,
                           tk::grm::statistics< Stack, use >,
                           tk::grm::pdfs< Stack, use, store_option > > > {};

  //! main keywords
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use >, walker > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore > {};

} // deck::
} // walker::

#endif // WalkerInputDeckGrammar_h
