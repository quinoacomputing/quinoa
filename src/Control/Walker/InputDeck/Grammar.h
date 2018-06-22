// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Walker's input deck grammar definition
  \details   Walker's input deck grammar definition. We use the Parsing
    Expression Grammar Template Library (PEGTL) to create the grammar and the
    associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef WalkerInputDeckGrammar_h
#define WalkerInputDeckGrammar_h

#include "Macro.h"
#include "Exception.h"
#include "Walker/Types.h"
#include "Keywords.h"
#include "CommonGrammar.h"
#include "QuinoaConfig.h"
#include "Walker/Options/InitPolicy.h"
#include "Walker/Options/CoeffPolicy.h"
#include "Walker/Options/HydroTimeScales.h"
#include "Walker/Options/HydroProductions.h"

#ifdef HAS_MKL
  #include "MKLGrammar.h"
#endif
#ifdef HAS_RNGSSE2
  #include "RNGSSEGrammar.h"
#endif
#include "Random123Grammar.h"

namespace walker {

extern ctr::InputDeck g_inputdeck_defaults;

//! Walker input deck facilitating user input for integrating SDEs
namespace deck {

  //! \brief Specialization of tk::grm::use for Walker's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword
                          , ctr::InputDeck::keywords1
                          , ctr::InputDeck::keywords2
                          , ctr::InputDeck::keywords3
                          , ctr::InputDeck::keywords4
                          , ctr::InputDeck::keywords5
                          , ctr::InputDeck::keywords6
                          , ctr::InputDeck::keywords7
                          , ctr::InputDeck::keywords8
                          >;

  // Walker's InputDeck state
 
  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  static tk::tuple::tagged_tuple< tag::dirichlet,       std::size_t,
                                  tag::gendir,          std::size_t,
                                  tag::wrightfisher,    std::size_t,
                                  tag::ou,              std::size_t,
                                  tag::diagou,          std::size_t,
                                  tag::skewnormal,      std::size_t,
                                  tag::gamma,           std::size_t,
                                  tag::velocity,        std::size_t,
                                  tag::position,        std::size_t,
                                  tag::dissipation,     std::size_t,
                                  tag::beta,            std::size_t,
                                  tag::numfracbeta,     std::size_t,
                                  tag::massfracbeta,    std::size_t,
                                  tag::mixnumfracbeta,  std::size_t,
                                  tag::mixmassfracbeta, std::size_t > neq;
} // ::deck
} // ::walker

namespace tk {
namespace grm {

  using namespace tao;

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // Walker's InputDeck actions

  //! Rule used to trigger action
  template< class Option, typename... tags >
  struct store_walker_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for walker.
  template< class Option, typename... tags >
  struct action< store_walker_option< Option, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      store_option< Stack, walker::deck::use, Option, walker::ctr::InputDeck,
                    Input, tags... >
                  ( stack, in, walker::g_inputdeck_defaults );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct register_eq : pegtl::success {};
  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  template< class eq >
  struct action< register_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {
      using walker::deck::neq;
      ++neq.get< eq >();
    }
  };

  //! Rule used to trigger action
  template< class eq, class vec, MsgKey key >
  struct check_vector_exists : pegtl::success {};
  //! \brief Check the existence of a vector (required for a block)
  //! \details This functor can be used to check the existence of a
  //!   user-specified vector, e.g., at the end of a diffeq ... end block, that
  //!   is required for that block. If the vector does not exist, we error out.
  //! \note This functor only checks existence. If the vector exists, the size
  //!   of it can be checked by check_vector_size.
  template< class eq, class vec, MsgKey Key >
  struct action< check_vector_exists< eq, vec, Key > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& vv = stack.template get< tag::param, eq, vec >();
      using walker::deck::neq;
      if (vv.size() != neq.get< eq >())
        Message< Stack, ERROR, Key >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq, class vec > struct check_vector_size : pegtl::success {};
  //! \brief Do error checking of a vector (required for a block)
  //! \details This functor can be used to verify the correct size of an already
  //!   existing vector, specified by the user in a given block. The vector is
  //!   required to be non-empty.
  //! \note This functor does not check existence of a vector. If the vector
  //!   does not even exist, it will throw an exception in DEBUG mode, while in
  //!   RELEASE mode it will attempt to access unallocated memory yielding a
  //!   segfault. The existence of the vector should be checked by
  //!   check_vector_exists first.
  template< class eq, class vec  >
  struct action< check_vector_size< eq, vec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& vv = stack.template get< tag::param, eq, vec >();
      Assert( !vv.empty(), "Vector of vectors checked must not be empty" );
      const auto& v = vv.back();
      if (v.empty())
        Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_eq : pegtl::success {};
  //! \brief Do general error checking on the differential equation block
  //! \details This is error checking that all equation types must satisfy.
  template< class eq >
  struct action< check_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      // Error out if no dependent variable has been selected
      const auto& depvar =
        stack.template get< tag::param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );

      // Error out if no number of components has been selected
      const auto& ncomp = stack.template get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NONCOMP >( stack, in );

      // Error out if no RNG has been selected
      const auto& rng = stack.template get< tag::param, eq, tag::rng >();
      if (rng.empty() || rng.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NORNG >( stack, in );

      // Error out if no initialization policy has been selected
      const auto& init =
        stack.template get< tag::param, eq, tag::initpolicy >();
      if (init.empty() || init.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOINIT >( stack, in );

      // Error out if no coefficients policy has been selected
      const auto& coeff =
        stack.template get< tag::param, eq, tag::coeffpolicy >();
      if (coeff.empty() || coeff.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOCOEFF >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_init : pegtl::success {};
  //! \brief Do error checking on the selected initialization policy
  template< class eq >
  struct action< check_init< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      const auto& init =
        stack.template get< tag::param, eq, tag::initpolicy >();
      // Error checks for joint delta initpolicy
      if (init.size() == neq.get< eq >() &&
          init.back() == walker::ctr::InitPolicyType::JOINTDELTA) {
        // Make sure there was an icdelta...end block with at least a single
        // spike...end block
        const auto& spike = stack.template get< tag::param, eq, tag::spike >();
        if (!spike.empty() && spike.back().empty())
          Message< Stack, ERROR, MsgKey::NODELTA >( stack, in );
      }
      // Error checks for joint beta initpolicy
      if (init.size() == neq.get< eq >() &&
          init.back() == walker::ctr::InitPolicyType::JOINTBETA) {
        // Make sure there was an icbeta...end block with at least a single
        // betapdf...end block
        const auto& betapdf =
          stack.template get< tag::param, eq, tag::betapdf >();
        if (!betapdf.empty() && betapdf.back().empty())
          Message< Stack, ERROR, MsgKey::NOBETA >( stack, in );
      }
      // Error checks for joint gamma initpolicy
      if (init.size() == neq.get< eq >() &&
          init.back() == walker::ctr::InitPolicyType::JOINTGAMMA) {
        // Make sure there was an icgamma...end block with at least a single
        // gammapdf...end block
        const auto& gammapdf =
          stack.template get< tag::param, eq, tag::gamma >();
        if (!gammapdf.empty() && gammapdf.back().empty())
          Message< Stack, ERROR, MsgKey::NOGAMMA >( stack, in );
      }
    }
  };

  //! Do error checking on eq coupled to another eq and compute depvar id
  //! \tparam eq Equation tag of equation
  //! \tparam coupledeq Equation tag of equation coupled to eq
  //! \tparam id Equation offsets tag for coupled equation
  //! \tparam depvar Error message key to use on missing coupled depvar
  //! \param[in] in Parser input
  //! \param[in,out] stack Grammar stack to wrok with
  //! \param[in] missing Error message key to use on missing coupled equation 
  template< typename eq, typename coupledeq, typename id, MsgKey depvar,
            typename Input, typename Stack >
  static void check_coupled( const Input& in, Stack& stack, MsgKey missing )
  {
    using walker::deck::neq;
    // Error out if a <eq> model is configured without a <coupledeq> model
    if (neq.get< coupledeq >() == 0) {
      stack.template push_back< tag::error >
        ( std::string("Parser error: ") + tk::cref_find( message, missing ) );
    } else {
      // get coupled position eq configuration
      const auto& ceq = stack.template get< tag::param, eq, coupledeq >();
      if (ceq.empty())
        // Error out if coupled <coupledeq> model eq depvar is not selected
        Message< Stack, ERROR, depvar >( stack, in );
      else { // find offset (local eq system index among systems) for depvar
        // get ncomponents object from this input deck
        const auto& ncomps = stack.template get< tag::component >();
        // compute offset map associating offsets to dependent variables
        auto offsetmap = ncomps.offsetmap( stack );
        // get and save offsets for all depvars for velocity eqs configured
        for (auto p : ceq)
          stack.template
            get< tag::param, eq, id >().
              push_back( tk::cref_find( offsetmap, p ) );
      }
    }
  }

  //! Rule used to trigger action
  struct check_velocity : pegtl::success {};
  //! \brief Do error checking on the velocity eq block
  template<>
  struct action< check_velocity > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Ensure a coupled position model is configured
      check_coupled< tag::velocity,
          tag::position, tag::position_id, MsgKey::POSITION_DEPVAR >
        ( in, stack, MsgKey::POSITION_MISSING );
      // Ensure a coupled dissipation model is configured
      check_coupled< tag::velocity,
          tag::dissipation, tag::dissipation_id, MsgKey::DISSIPATION_DEPVAR >
        ( in, stack, MsgKey::DISSIPATION_MISSING );
    }
  };

  //! Rule used to trigger action
  struct check_position : pegtl::success {};
  //! \brief Do error checking on the position eq block
  template<>
  struct action< check_position > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Ensure a coupled velocity model is configured
      check_coupled< tag::position,
          tag::velocity, tag::velocity_id, MsgKey::VELOCITY_DEPVAR >
        ( in, stack, MsgKey::VELOCITY_MISSING );
    }
  };

  //! Rule used to trigger action
  struct check_dissipation : pegtl::success {};
  //! \brief Do error checking on the dissipation eq block
  template<>
  struct action< check_dissipation > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Ensure a coupled velocity model is configured
      check_coupled< tag::dissipation,
          tag::velocity, tag::velocity_id, MsgKey::VELOCITY_DEPVAR >
        ( in, stack, MsgKey::VELOCITY_MISSING );
    }
  };

  //! Rule used to trigger action
  struct velocity_defaults : pegtl::success {};
  //! \brief Set defaults for the velocity model
  template<>
  struct action< velocity_defaults > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      using walker::deck::neq;
      // Set number of components: always 3 velocity components
      auto& ncomp = stack.template get< tag::component, tag::velocity >();
      ncomp.push_back( 3 );
      // Set C0 = 2.1 if not specified
      auto& C0 = stack.template get< tag::param, tag::velocity, tag::c0 >();
      if (C0.size() != neq.get< tag::velocity >()) C0.push_back( 2.1 );
    }
  };

  //! Rule used to trigger action
  struct position_defaults : pegtl::success {};
  //! \brief Set defaults for all Lagrangian particle position models
  template<>
  struct action< position_defaults > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      using walker::deck::neq;
      // Set number of components: always 3 position components
      auto& ncomp = stack.template get< tag::component, tag::position >();
      ncomp.push_back( 3 );
      // Set RNG if no RNG has been selected (not all position models use this,
      // so don't impose on the user to define one). Pick a Random123 generator,
      // as that is always available.
      auto& rngs = stack.template get< tag::selected, tag::rng >();
      auto& rng = stack.template get< tag::param, tag::position, tag::rng >();
      if (rng.empty() || rng.size() != neq.get< tag::position >()) {
        // add RNG to the list of selected RNGs (as it has been parsed)
        rngs.push_back( tk::ctr::RNG().value( kw::r123_philox::string() ) );
        // select RNG for position equation (as it has been parsed)
        rng.push_back( tk::ctr::RNGType::R123_PHILOX );
      }
    }
  };

  //! Rule used to trigger action
  struct dissipation_defaults : pegtl::success {};
  //! \brief Set defaults for all Lagrangian particle dissipation models
  template<>
  struct action< dissipation_defaults > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      using walker::deck::neq;
      // Set number of components: always 1 dissipation component
      auto& ncomp = stack.template get< tag::component, tag::dissipation >();
      ncomp.push_back( 1 );
    }
  };

} // ::grm
} // ::tk

namespace walker {

//! Walker input deck facilitating user input for integrating SDEs
namespace deck {

  using namespace tao;

  // Walker's InputDeck grammar

  //! scan and store_back sde keyword and option
  template< typename keyword, class eq >
  struct scan_sde :
         tk::grm::scan< typename keyword::pegtl_string,
                        tk::grm::store_back_option< use,
                                                    ctr::DiffEq,
                                                    tag::selected,
                                                    tag::diffeq >,
                        // start new vector or vectors of spikes for a potential
                        // jointdelta initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::spike >,
                        // start new vector or vectors of beta parameters for a
                        // potential jointbeta initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::betapdf >,
                        // start new vector or vectors of gamma parameters for a
                        // potential jointgamma initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::gamma >,
                        // start new vector or vectors of gaussian parameters
                        // for a potential jointgaussian initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::gaussian > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discrparam< use, kw::npar, tag::npar >,
                     tk::grm::discrparam< use, kw::nstep, tag::nstep >,
                     tk::grm::discrparam< use, kw::term, tag::term >,
                     tk::grm::discrparam< use, kw::dt, tag::dt >,
                     tk::grm::interval< use< kw::ttyi >, tag::tty > > {};

  //! rngs
  struct rngs :
         pegtl::sor<
                     #ifdef HAS_MKL
                     tk::mkl::rngs< use,
                                    tag::selected, tag::rng,
                                    tag::param, tag::rngmkl >,
                     #endif
                     #ifdef HAS_RNGSSE2
                     tk::rngsse::rngs< use,
                                       tag::selected, tag::rng,
                                       tag::param, tag::rngsse >,
                     #endif
                     tk::random123::rngs< use,
                                          tag::selected, tag::rng,
                                          tag::param, tag::rng123 > > {};

  //! scan icdelta ... end block
  template< class eq >
  struct icdelta :
         pegtl::if_must<
           tk::grm::readkw< use< kw::icdelta >::pegtl_string >,
           // parse a spike ... end block (there can be multiple)
           tk::grm::block< use< kw::end >,
                           tk::grm::parameter_vector<
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
         pegtl::if_must<
           tk::grm::readkw< use< kw::icbeta >::pegtl_string >,
           // parse a betapdf ... end block (there can be multiple)
           tk::grm::block< use< kw::end >,
                           tk::grm::parameter_vector<
                             use,
                             use< kw::betapdf >,
                             tk::grm::Store_back_back_back,
                             tk::grm::start_vector_back,
                             tk::grm::check_betapdfs,
                             eq,
                             tag::betapdf > > > {};

  //! scan icgamma ... end block
  template< class eq >
  struct icgamma :
         pegtl::if_must<
           tk::grm::readkw< use< kw::icgamma >::pegtl_string >,
           // parse a gammapdf ... end block (there can be multiple)
           tk::grm::block< use< kw::end >,
                           tk::grm::parameter_vector<
                             use,
                             use< kw::gammapdf >,
                             tk::grm::Store_back_back_back,
                             tk::grm::start_vector_back,
                             tk::grm::check_gammapdfs,
                             eq,
                             tag::gamma > > > {};

  //! scan icgaussian ... end block
  template< class eq >
  struct icgaussian :
         pegtl::if_must<
           tk::grm::readkw< use< kw::icgaussian >::pegtl_string >,
           // parse a gaussian ... end block (there can be multiple)
           tk::grm::block< use< kw::end >,
                           tk::grm::parameter_vector<
                             use,
                             use< kw::gaussian >,
                             tk::grm::Store_back_back_back,
                             tk::grm::start_vector_back,
                             tk::grm::check_gaussians,
                             eq,
                             tag::gaussian > > > {};

  //! Error checks after an equation ... end block has been parsed
  template< class eq, class... extra_checks  >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           tk::grm::register_eq< eq >,
           // do error checking on this block
           tk::grm::check_eq< eq >,
           // do error checking on the init policy
           tk::grm::check_init< eq >,
           // performe extra pegtl actions, e.g., performing extra checks
           extra_checks... > {};

  //! SDE parameter vector
  template< class keyword, class eq, class param >
  struct sde_parameter_vector :
         tk::grm::parameter_vector< use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq,
                                    param > {};

  //! SDE option vector
  template< class Option, class keyword, class eq, class param,
            template< class, class > class check >
  struct sde_option_vector :
         tk::grm::option_vector< use,
                                 use< keyword >,
                                 Option,
                                 tk::grm::Store_back_back,
                                 tk::grm::start_vector,
                                 check,
                                 eq,
                                 param > {};

  //! Diagonal Ornstein-Uhlenbeck SDE
  struct diag_ou :
         pegtl::if_must<
           scan_sde< use< kw::diag_ou >, tag::diagou >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::diagou,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::diagou >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::diagou,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::diagou,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::diagou,
                                            tag::coeffpolicy >,
                           icdelta< tag::diagou >,
                           icbeta< tag::diagou >,
                           icgamma< tag::diagou >,
                           icgaussian< tag::diagou >,
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
         pegtl::if_must<
           scan_sde< use< kw::ornstein_uhlenbeck >, tag::ou >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::ou,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::ou >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::ou,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::ou,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::ou,
                                            tag::coeffpolicy >,
                           icdelta< tag::ou >,
                           icbeta< tag::ou >,
                           icgamma< tag::ou >,
                           icgaussian< tag::ou >,
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
         pegtl::if_must<
           scan_sde< use< kw::skewnormal >, tag::skewnormal >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::skewnormal,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::skewnormal >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::skewnormal,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::skewnormal,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::skewnormal,
                                            tag::coeffpolicy >,
                           icdelta< tag::skewnormal >,
                           icbeta< tag::skewnormal >,
                           icgamma< tag::skewnormal >,
                           icgaussian< tag::skewnormal >,
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
         pegtl::if_must<
           scan_sde< use< kw::beta >, tag::beta >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::beta,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::beta >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::beta,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::beta,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::beta,
                                            tag::coeffpolicy >,
                           icdelta< tag::beta >,
                           icbeta< tag::beta >,
                           icgamma< tag::beta >,
                           icgaussian< tag::beta >,
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
         pegtl::if_must<
           scan_sde< use< kw::numfracbeta >, tag::numfracbeta >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::numfracbeta,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::numfracbeta >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::numfracbeta,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::numfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::numfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::numfracbeta >,
                           icbeta< tag::numfracbeta >,
                           icgamma< tag::numfracbeta >,
                           icgaussian< tag::numfracbeta >,
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
         pegtl::if_must<
           scan_sde< use< kw::massfracbeta >, tag::massfracbeta >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::massfracbeta,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::massfracbeta >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::massfracbeta,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::massfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::massfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::massfracbeta >,
                           icbeta< tag::massfracbeta >,
                           icgamma< tag::massfracbeta >,
                           icgaussian< tag::massfracbeta >,
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
         pegtl::if_must<
           scan_sde< use< kw::mixnumfracbeta >, tag::mixnumfracbeta >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::mixnumfracbeta,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::mixnumfracbeta >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::mixnumfracbeta,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::mixnumfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::mixnumfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::mixnumfracbeta >,
                           icbeta< tag::mixnumfracbeta >,
                           icgamma< tag::mixnumfracbeta >,
                           icgaussian< tag::mixnumfracbeta >,
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
         pegtl::if_must<
           scan_sde< use< kw::mixmassfracbeta >, tag::mixmassfracbeta >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::mixmassfracbeta,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::mixmassfracbeta >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::mixmassfracbeta,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::mixmassfracbeta,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::mixmassfracbeta,
                                            tag::coeffpolicy >,
                           icdelta< tag::mixmassfracbeta >,
                           icbeta< tag::mixmassfracbeta >,
                           icgamma< tag::mixmassfracbeta >,
                           icgaussian< tag::mixmassfracbeta >,
                           sde_option_vector< ctr::HydroTimeScales,
                                              kw::hydrotimescales,
                                              tag::mixmassfracbeta,
                                              tag::hydrotimescales,
                                              tk::grm::check_vector_size >,
                           sde_option_vector< ctr::HydroProductions,
                                              kw::hydroproductions,
                                              tag::mixmassfracbeta,
                                              tag::hydroproductions,
                                              tk::grm::check_vector_size >,
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
           check_errors< tag::mixmassfracbeta,
                         tk::grm::check_vector_exists<
                           tag::mixmassfracbeta,
                           tag::hydrotimescales,
                           tk::grm::MsgKey::HYDROTIMESCALES >,
                         tk::grm::check_vector_exists<
                           tag::mixmassfracbeta,
                           tag::hydroproductions,
                           tk::grm::MsgKey::HYDROPRODUCTIONS > > > {};

  //! Gamma SDE
  struct gamma :
         pegtl::if_must<
           scan_sde< use< kw::gamma >, tag::gamma >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::gamma,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::gamma >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::gamma,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::gamma,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::gamma,
                                            tag::coeffpolicy >,
                           icdelta< tag::gamma >,
                           icbeta< tag::gamma >,
                           icgamma< tag::gamma >,
                           icgaussian< tag::gamma >,
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
         pegtl::if_must<
           scan_sde< use< kw::dirichlet >, tag::dirichlet >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::dirichlet,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::dirichlet >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::dirichlet,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::dirichlet,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::dirichlet,
                                            tag::coeffpolicy >,
                           icdelta< tag::dirichlet >,
                           icbeta< tag::dirichlet >,
                           icgamma< tag::dirichlet >,
                           icgaussian< tag::dirichlet >,
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
         pegtl::if_must<
           scan_sde< use< kw::gendir >, tag::gendir >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::gendir,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::gendir >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::gendir,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::gendir,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::gendir,
                                            tag::coeffpolicy >,
                           icdelta< tag::gendir >,
                           icbeta< tag::gendir >,
                           icgamma< tag::gendir >,
                           icgaussian< tag::gendir >,
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
         pegtl::if_must<
           scan_sde< use< kw::wrightfisher >, tag::wrightfisher >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::wrightfisher,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::wrightfisher >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::wrightfisher,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::wrightfisher,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::wrightfisher,
                                            tag::coeffpolicy >,
                           icdelta< tag::wrightfisher >,
                           icbeta< tag::wrightfisher >,
                           icgamma< tag::wrightfisher >,
                           icgaussian< tag::wrightfisher >,
                           sde_parameter_vector< kw::sde_omega,
                                                 tag::wrightfisher,
                                                 tag::omega > >,
           check_errors< tag::wrightfisher > > {};

  //! Velocity SDE
  struct velocity :
         pegtl::if_must<
           scan_sde< use< kw::velocity >, tag::velocity >,
           tk::grm::velocity_defaults,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::velocity,
                                            tag::depvar >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::velocity,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::velocity,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::velocity,
                                            tag::coeffpolicy >,
                           icdelta< tag::velocity >,
                           icbeta< tag::velocity >,
                           icgamma< tag::velocity >,
                           icgaussian< tag::velocity >,
                           sde_option_vector< ctr::HydroTimeScales,
                                              kw::hydrotimescales,
                                              tag::velocity,
                                              tag::hydrotimescales,
                                              tk::grm::check_vector_size >,
                           sde_option_vector< ctr::HydroProductions,
                                              kw::hydroproductions,
                                              tag::velocity,
                                              tag::hydroproductions,
                                              tk::grm::check_vector_size >,
                           tk::grm::process<
                             use< kw::sde_c0 >,
                             tk::grm::Store_back< tag::param,
                                                  tag::velocity,
                                                  tag::c0 > >,
                           tk::grm::process<
                             use< kw::position >,
                             tk::grm::Store_back< tag::param,
                                                  tag::velocity,
                                                  tag::position >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::dissipation >,
                             tk::grm::Store_back< tag::param,
                                                  tag::velocity,
                                                  tag::dissipation >,
                             pegtl::alpha > >,
           check_errors< tag::velocity > > {};

  //! position equation
  struct position :
         pegtl::if_must<
           scan_sde< use< kw::position >, tag::position >,
           tk::grm::position_defaults,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::position,
                                            tag::depvar >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::position,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::position,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::position,
                                            tag::coeffpolicy >,
                           icdelta< tag::position >,
                           icbeta< tag::position >,
                           icgamma< tag::position >,
                           icgaussian< tag::position >,
                           tk::grm::process<
                             use< kw::velocity >,
                             tk::grm::Store_back< tag::param,
                                                  tag::position,
                                                  tag::velocity >,
                             pegtl::alpha > >,
           check_errors< tag::position > > {};

  //! dissipation equation
  struct dissipation :
         pegtl::if_must<
           scan_sde< use< kw::dissipation >, tag::dissipation >,
           tk::grm::dissipation_defaults,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::dissipation,
                                            tag::depvar >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::dissipation,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::dissipation,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::dissipation,
                                            tag::coeffpolicy >,
                           icdelta< tag::dissipation >,
                           icbeta< tag::dissipation >,
                           icgamma< tag::dissipation >,
                           icgaussian< tag::dissipation >,
                           tk::grm::process<
                             use< kw::velocity >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::velocity >,
                             pegtl::alpha > >,
           check_errors< tag::dissipation > > {};

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
                     mixmassfracbeta,
                     position,
                     dissipation,
                     velocity > {};


  //! 'walker' block
  struct walker :
         pegtl::if_must<
           tk::grm::readkw< use< kw::walker >::pegtl_string >,
           pegtl::sor<
             pegtl::seq<
               tk::grm::block<
                 use< kw::end >,
                 discretization_parameters,
                 sde,
                 tk::grm::rngblock< use, rngs >,
                 tk::grm::statistics< use, tk::grm::store_walker_option >,
                 tk::grm::pdfs< use, tk::grm::store_walker_option > >,
               tk::grm::check_velocity,
               tk::grm::check_position,
               tk::grm::check_dissipation >,
           tk::grm::msg< tk::grm::MsgType::ERROR,
                         tk::grm::MsgKey::UNFINISHED > > > {};

  //! main keywords
  struct keywords :
         pegtl::sor< tk::grm::title< use >, walker > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< keywords, tk::grm::ignore > {};

} // deck::
} // walker::

#endif // WalkerInputDeckGrammar_h
