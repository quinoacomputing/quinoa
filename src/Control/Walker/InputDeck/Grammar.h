// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Tue 26 Jul 2016 07:40:52 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker's input deck grammar definition
  \details   Walker's input deck grammar definition. We use the [Parsing
    Expression Grammar Template Library (PEGTL)]
    (https://code.google.com/p/pegtl/wiki/PEGTL0) to create the grammar and the
    associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL.
    Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef WalkerInputDeckGrammar_h
#define WalkerInputDeckGrammar_h

#include "Macro.h"
#include "Exception.h"
#include "PEGTLParsed.h"
#include "Walker/Types.h"
#include "Keywords.h"
#include "CommonGrammar.h"
#include "QuinoaConfig.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/CoeffPolicy.h"

#ifdef HAS_MKL
#include "MKLGrammar.h"
#endif

#include "RNGSSEGrammar.h"

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
                            ctr::InputDeck::keywords5,
                            ctr::InputDeck::keywords6 >;

  // Walker's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  //! \author J. Bakosi
  static tk::tuple::tagged_tuple< tag::dirichlet,       std::size_t,
                                  tag::gendir,          std::size_t,
                                  tag::wrightfisher,    std::size_t,
                                  tag::ou,              std::size_t,
                                  tag::diagou,          std::size_t,
                                  tag::skewnormal,      std::size_t,
                                  tag::gamma,           std::size_t,
                                  tag::beta,            std::size_t,
                                  tag::numfracbeta,     std::size_t,
                                  tag::massfracbeta,    std::size_t,
                                  tag::mixnumfracbeta,  std::size_t,
                                  tag::mixmassfracbeta, std::size_t > neq;

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

  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  //! \author J. Bakosi
  template< class eq >
  struct register_eq : pegtl::action_base< register_eq< eq > > {
    static void apply( const std::string&, Stack& ) {
      ++neq.get< eq >();
    }
  };

  //! \brief Do general error checking on the differential equation block
  //! \details This is error checking that all equation types must satisfy.
  //! \author J. Bakosi
  template< class eq >
  struct check_eq : pegtl::action_base< check_eq< eq > > {
    static void apply( const std::string& value, Stack& stack ) {

      // Error out if no dependent variable has been selected
      const auto& depvar = stack.get< tag::param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NODEPVAR >
                        ( stack, value );

      // Error out if no number of components has been selected
      const auto& ncomp = stack.get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NONCOMP >
                        ( stack, value );

      // Error out if no RNG has been selected
      const auto& rng = stack.get< tag::param, eq, tag::rng >();
      if (rng.empty() || rng.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NORNG >
                        ( stack, value );

      // Error out if no initialization policy has been selected
      const auto& init = stack.get< tag::param, eq, tag::initpolicy >();
      if (init.empty() || init.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOINIT >
                        ( stack, value );

      // Error out if no coefficients policy has been selected
      const auto& coeff = stack.get< tag::param, eq, tag::coeffpolicy >();
      if (coeff.empty() || coeff.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOCOEFF >
                        ( stack, value );
    }
  };

  //! \brief Do error checking on the selected initialization policy
  //! \author J. Bakosi
  template< class eq >
  struct check_init : pegtl::action_base< check_init< eq > > {
    static void apply( const std::string& value, Stack& stack ) {
      const auto& init = stack.get< tag::param, eq, tag::initpolicy >();
      // Error checks for joint delta initpolicy
      if (init.size() == neq.get< eq >() &&
          init.back() == ctr::InitPolicyType::JOINTDELTA) {
        // Make sure there was an icdelta...end block with at least a single
        // spike...end block
        const auto& spike = stack.template get< tag::param, eq, tag::spike >();
        if (!spike.empty() && spike.back().empty())
          tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NODELTA >
                          ( stack, value );
      }
      // Error checks for joint beta initpolicy
      if (init.size() == neq.get< eq >() &&
          init.back() == ctr::InitPolicyType::JOINTBETA) {
        // Make sure there was an icbeta...end block with at least a single
        // betapdf...end block
        const auto& betapdf =
          stack.template get< tag::param, eq, tag::betapdf >();
        if (!betapdf.empty() && betapdf.back().empty())
          tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOBETA >
                          ( stack, value );
      }

    }
  };


  // Walker's InputDeck grammar

  //! scan and store_back sde keyword and option
  template< typename keyword, class eq >
  struct scan_sde :
         tk::grm::scan< Stack,
                        typename keyword::pegtl_string,
                        tk::grm::store_back_option< Stack,
                                                    use,
                                                    ctr::DiffEq,
                                                    tag::selected,
                                                    tag::diffeq >,
                        // start new vector or vectors of spikes for a potential
                        // jointdelta initpolicy
                        tk::grm::start_vector< Stack,
                                               tag::param,
                                               eq,
                                               tag::spike >,
                        // start new vector or vectors of beta parameters for a
                        // potential jointbeta initpolicy
                        tk::grm::start_vector< Stack,
                                               tag::param,
                                               eq,
                                               tag::betapdf > > {};

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

  //! scan icdelta ... end block
  template< class eq >
  struct icdelta :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::icdelta >::pegtl_string >,
           // parse a spike ... end block (there can be multiple)
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::parameter_vector<
                             Stack,
                             use,
                             use< kw::spike >,
                             tk::grm::Store_back_back_back,
                             tk::grm::start_vector_back,
                             tk::grm::check_spikes,
                             eq,
                             tag::spike > > > {};

  //! scan icbeta ... end block
  template< class eq >
  struct icbeta :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::icbeta >::pegtl_string >,
           // parse a betapdf ... end block (there can be multiple)
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::parameter_vector<
                             Stack,
                             use,
                             use< kw::betapdf >,
                             tk::grm::Store_back_back_back,
                             tk::grm::start_vector_back,
                             tk::grm::check_betapdfs,
                             eq,
                             tag::betapdf > > > {};

  //! Error checks after a equation..end block has been parsed
  template< class eq >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           pegtl::apply< register_eq< eq > >,
           // do error checking on this block
           pegtl::apply< check_eq< eq > >,
           // do error checking on the init policy
           pegtl::apply< check_init< eq > > > {};

  //! SDE parameter vector
  template< class keyword, class eq, class param >
  struct sde_parameter_vector :
         tk::grm::parameter_vector< Stack,
                                    use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq,
                                    param > {};

  //! Diagonal Ornstein-Uhlenbeck SDE
  struct diag_ou :
         pegtl::ifmust<
           scan_sde< use< kw::diag_ou >, tag::diagou >,
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
                                            ctr::InitPolicy,
                                            tag::diagou,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::diagou,
                                            tag::coeffpolicy >,
                           icdelta< tag::diagou >,
                           icbeta< tag::diagou >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tag::diagou,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_theta,
                                                 tag::diagou,
                                                 tag::theta >,
                           sde_parameter_vector< kw::sde_mu,
                                                 tag::diagou,
                                                 tag::mu > >,
           check_errors< tag::diagou > > {};

  //! Ornstein-Uhlenbeck SDE
  struct ornstein_uhlenbeck :
         pegtl::ifmust<
           scan_sde< use< kw::ornstein_uhlenbeck >, tag::ou >,
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
                                            ctr::InitPolicy,
                                            tag::ou,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::ou,
                                            tag::coeffpolicy >,
                           icdelta< tag::ou >,
                           icbeta< tag::ou >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tag::ou,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_theta,
                                                 tag::ou,
                                                 tag::theta >,
                           sde_parameter_vector< kw::sde_mu,
                                                 tag::ou,
                                                 tag::mu > >,
           check_errors< tag::ou > > {};

  //! Skew-normal SDE
  struct skewnormal :
         pegtl::ifmust<
           scan_sde< use< kw::skewnormal >, tag::skewnormal >,
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
                                            ctr::InitPolicy,
                                            tag::skewnormal,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::skewnormal,
                                            tag::coeffpolicy >,
                           icdelta< tag::skewnormal >,
                           icbeta< tag::skewnormal >,
                           sde_parameter_vector< kw::sde_T,
                                                 tag::skewnormal,
                                                 tag::timescale >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tag::skewnormal,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_lambda,
                                                 tag::skewnormal,
                                                 tag::lambda > >,
           check_errors< tag::skewnormal > > {};

  //! Beta SDE
  struct beta :
         pegtl::ifmust<
           scan_sde< use< kw::beta >, tag::beta >,
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
                                            ctr::InitPolicy,
                                            tag::beta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::beta,
                                            tag::coeffpolicy >,
                           icdelta< tag::beta >,
                           icbeta< tag::beta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::beta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::beta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::beta,
                                                 tag::kappa > >,
           check_errors< tag::beta > > {};

  //! Number-fraction beta SDE
  struct numfracbeta :
         pegtl::ifmust<
           scan_sde< use< kw::numfracbeta >, tag::numfracbeta >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::depvar< Stack,
                                            use,
                                            tag::numfracbeta,
                                            tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::numfracbeta >,
                           tk::grm::rng< Stack,
                                         use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::numfracbeta,
                                         tag::rng >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::numfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::numfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::numfracbeta >,
                           icbeta< tag::numfracbeta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::numfracbeta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::numfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::numfracbeta,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tag::numfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_rcomma,
                                                 tag::numfracbeta,
                                                 tag::rcomma > >,
           check_errors< tag::numfracbeta > > {};

  //! Mass-fraction beta SDE
  struct massfracbeta :
         pegtl::ifmust<
           scan_sde< use< kw::massfracbeta >, tag::massfracbeta >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::depvar< Stack,
                                            use,
                                            tag::massfracbeta,
                                            tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::massfracbeta >,
                           tk::grm::rng< Stack,
                                         use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::massfracbeta,
                                         tag::rng >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::massfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::massfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::massfracbeta >,
                           icbeta< tag::massfracbeta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::massfracbeta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::massfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::massfracbeta,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tag::massfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_r,
                                                 tag::massfracbeta,
                                                 tag::r > >,
           check_errors< tag::massfracbeta > > {};

  //! Mix number-fraction beta SDE
  struct mixnumfracbeta :
         pegtl::ifmust<
           scan_sde< use< kw::mixnumfracbeta >, tag::mixnumfracbeta >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::depvar< Stack,
                                            use,
                                            tag::mixnumfracbeta,
                                            tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::mixnumfracbeta >,
                           tk::grm::rng< Stack,
                                         use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::mixnumfracbeta,
                                         tag::rng >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::mixnumfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::mixnumfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::mixnumfracbeta >,
                           icbeta< tag::mixnumfracbeta >,
                           sde_parameter_vector< kw::sde_bprime,
                                                 tag::mixnumfracbeta,
                                                 tag::bprime >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::mixnumfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappaprime,
                                                 tag::mixnumfracbeta,
                                                 tag::kappaprime >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tag::mixnumfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_rcomma,
                                                 tag::mixnumfracbeta,
                                                 tag::rcomma > >,
           check_errors< tag::mixnumfracbeta > > {};

  //! Mix mass-fraction beta SDE
  struct mixmassfracbeta :
         pegtl::ifmust<
           scan_sde< use< kw::mixmassfracbeta >, tag::mixmassfracbeta >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::depvar< Stack,
                                            use,
                                            tag::mixmassfracbeta,
                                            tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::mixmassfracbeta >,
                           tk::grm::rng< Stack,
                                         use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::mixmassfracbeta,
                                         tag::rng >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::mixmassfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::mixmassfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::mixmassfracbeta >,
                           icbeta< tag::mixmassfracbeta >,
                           sde_parameter_vector< kw::sde_bprime,
                                                 tag::mixmassfracbeta,
                                                 tag::bprime >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::mixmassfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappaprime,
                                                 tag::mixmassfracbeta,
                                                 tag::kappaprime >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tag::mixmassfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_r,
                                                 tag::mixmassfracbeta,
                                                 tag::r > >,
           check_errors< tag::mixmassfracbeta > > {};

  //! Gamma SDE
  struct gamma :
         pegtl::ifmust<
           scan_sde< use< kw::gamma >, tag::gamma >,
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
                                            ctr::InitPolicy,
                                            tag::gamma,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::gamma,
                                            tag::coeffpolicy >,
                           icdelta< tag::gamma >,
                           icbeta< tag::gamma >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::gamma,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::gamma,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::gamma,
                                                 tag::kappa > >,
           check_errors< tag::gamma > > {};

  //! Dirichlet SDE
  struct dirichlet :
         pegtl::ifmust<
           scan_sde< use< kw::dirichlet >, tag::dirichlet >,
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
                                            ctr::InitPolicy,
                                            tag::dirichlet,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::dirichlet,
                                            tag::coeffpolicy >,
                           icdelta< tag::dirichlet >,
                           icbeta< tag::dirichlet >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::dirichlet,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::dirichlet,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::dirichlet,
                                                 tag::kappa > >,
           check_errors< tag::dirichlet > > {};

  //! Generalized Dirichlet SDE
  struct gendir :
         pegtl::ifmust<
           scan_sde< use< kw::gendir >, tag::gendir >,
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
                                            ctr::InitPolicy,
                                            tag::gendir,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::gendir,
                                            tag::coeffpolicy >,
                           icdelta< tag::gendir >,
                           icbeta< tag::gendir >,
                           sde_parameter_vector< kw::sde_b,
                                                 tag::gendir,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tag::gendir,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tag::gendir,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_c,
                                                 tag::gendir,
                                                 tag::c > >,
           check_errors< tag::gendir > > {};

  //! Wright-Fisher SDE
  struct wright_fisher :
         pegtl::ifmust<
           scan_sde< use< kw::wrightfisher >, tag::wrightfisher >,
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
                                            ctr::InitPolicy,
                                            tag::wrightfisher,
                                            tag::initpolicy >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::wrightfisher,
                                            tag::coeffpolicy >,
                           icdelta< tag::wrightfisher >,
                           icbeta< tag::wrightfisher >,
                           sde_parameter_vector< kw::sde_omega,
                                                 tag::wrightfisher,
                                                 tag::omega > >,
           check_errors< tag::wrightfisher > > {};

  //! stochastic differential equations
  struct sde :
         pegtl::sor< dirichlet,
                     gendir,
                     wright_fisher,
                     ornstein_uhlenbeck,
                     diag_ou,
                     skewnormal,
                     gamma,
                     beta,
                     numfracbeta,
                     massfracbeta,
                     mixnumfracbeta,
                     mixmassfracbeta > {};

  //! 'walker' block
  struct walker :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::walker >::pegtl_string >,
           pegtl::sor< tk::grm::block< Stack,
                         use< kw::end >,
                         discretization_parameters,
                         sde,
                         tk::grm::rngblock< Stack, use, rngs >,
                         tk::grm::statistics< Stack, use, store_option >,
                         tk::grm::pdfs< Stack, use, store_option > >,
                       pegtl::apply<
                          tk::grm::error< Stack,
                                          tk::grm::MsgKey::UNFINISHED > > > > {};

  //! main keywords
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use >, walker > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore< Stack > > {};

} // deck::
} // walker::

#endif // WalkerInputDeckGrammar_h
