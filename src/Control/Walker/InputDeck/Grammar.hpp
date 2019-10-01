// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Walker's input deck grammar definition
  \details   Walker's input deck grammar definition. We use the Parsing
    Expression Grammar Template Library (PEGTL) to create the grammar and the
    associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef WalkerInputDeckGrammar_h
#define WalkerInputDeckGrammar_h

#include <limits>
#include <algorithm>

#include "Macro.hpp"
#include "Exception.hpp"
#include "Walker/Types.hpp"
#include "Keywords.hpp"
#include "CommonGrammar.hpp"
#include "QuinoaConfig.hpp"
#include "Walker/Options/InitPolicy.hpp"
#include "Walker/Options/CoeffPolicy.hpp"
#include "Walker/Options/HydroTimeScales.hpp"
#include "Walker/Options/HydroProductions.hpp"

#ifdef HAS_MKL
  #include "MKLGrammar.hpp"
#endif
#ifdef HAS_RNGSSE2
  #include "RNGSSEGrammar.hpp"
#endif
#include "Random123Grammar.hpp"

namespace walker {

extern ctr::InputDeck g_inputdeck_defaults;

//! Walker input deck facilitating user input for integrating SDEs
namespace deck {

  //! \brief Specialization of tk::grm::use for Walker's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::InputDeck::keywords >;

  // Walker's InputDeck state
 
  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  static tk::TaggedTuple< brigand::list<
      tag::dirichlet,       std::size_t
    , tag::mixdirichlet,    std::size_t
    , tag::gendir,          std::size_t
    , tag::wrightfisher,    std::size_t
    , tag::ou,              std::size_t
    , tag::diagou,          std::size_t
    , tag::skewnormal,      std::size_t
    , tag::gamma,           std::size_t
    , tag::velocity,        std::size_t
    , tag::position,        std::size_t
    , tag::dissipation,     std::size_t
    , tag::beta,            std::size_t
    , tag::numfracbeta,     std::size_t
    , tag::massfracbeta,    std::size_t
    , tag::mixnumfracbeta,  std::size_t
    , tag::mixmassfracbeta, std::size_t
  > > neq;

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
  template< class eq, class vec >
  struct check_mean_gradient : pegtl::success {};
  //! Do error checking for a vector of prescribed mean gradient
  template< class eq, class vec  >
  struct action< check_mean_gradient< eq, vec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto& vv = stack.template get< tag::param, eq, vec >();
      Assert( !vv.empty(), "Vector of vectors checked must not be empty" );
      if (vv.back().size() != 3)
        Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq, class vec >
  struct check_gravity : pegtl::success {};
  //! Do error checking for a vector of prescribed mean gradient
  template< class eq, class vec  >
  struct action< check_gravity< eq, vec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      auto& vv = stack.template get< tag::param, eq, vec >();
      Assert( !vv.empty(), "Vector of vectors checked must not be empty" );
      if (vv.back().size() != 3)
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
      const auto& initpolicy =
        stack.template get< tag::param, eq, tag::initpolicy >();
      // Error checks for joint delta initpolicy
      if (initpolicy.size() == neq.get< eq >() &&
          initpolicy.back() == walker::ctr::InitPolicyType::JOINTDELTA) {
        // Make sure there was an icdelta...end block with at least a single
        // spike...end block
        const auto& spike =
          stack.template get< tag::param, eq, tag::init, tag::spike >();
        if (!spike.empty() && spike.back().empty())
          Message< Stack, ERROR, MsgKey::NODELTA >( stack, in );
      }
      // Error checks for joint beta initpolicy
      if (initpolicy.size() == neq.get< eq >() &&
          initpolicy.back() == walker::ctr::InitPolicyType::JOINTBETA) {
        // Make sure there was an icbeta...end block with at least a single
        // betapdf...end block
        const auto& betapdf =
          stack.template get< tag::param, eq, tag::init, tag::betapdf >();
        if (!betapdf.empty() && betapdf.back().empty())
          Message< Stack, ERROR, MsgKey::NOBETA >( stack, in );
      }
      // Error checks for joint gamma initpolicy
      if (initpolicy.size() == neq.get< eq >() &&
          initpolicy.back() == walker::ctr::InitPolicyType::JOINTGAMMA) {
        // Make sure there was an icgamma...end block with at least a single
        // gammapdf...end block
        const auto& gammapdf =
          stack.template get< tag::param, eq, tag::init, tag::gamma >();
        if (!gammapdf.empty() && gammapdf.back().empty())
          Message< Stack, ERROR, MsgKey::NOGAMMA >( stack, in );
      }
      // Error checks for joint correlated Gaussian initpolicy
      if (initpolicy.size() == neq.get< eq >() &&
          initpolicy.back() == walker::ctr::InitPolicyType::JOINTCORRGAUSSIAN) {
        // Ensure there was a mean vector and covaraiance matrix configured
        const auto& mean =
          stack.template get< tag::param, eq, tag::init, tag::mean >();
        if (mean.empty() || mean.back().empty())
          Message< Stack, ERROR, MsgKey::NOMEAN >( stack, in );
        const auto& cov =
          stack.template get< tag::param, eq, tag::init, tag::cov >();
        if (cov.empty() || cov.back().empty())
          Message< Stack, ERROR, MsgKey::NOCOV >( stack, in );
        // Ensure that an MKL RNG is configured if initpolicy is corr-Gaussian
        const auto& rng = stack.template get< tag::param, eq, tag::rng >();
        if (tk::ctr::RNG().lib( rng.back() ) != tk::ctr::RNGLibType::MKL)
          Message< Stack, ERROR, MsgKey::NOMKLRNG >( stack, in );
      }
    }
  };

  //! Setup coupling between two equations
  //! \tparam eq Tag of equation to be coupled
  //! \tparam coupledeq Tag of equation coupled to eq
  //! \tparam id Tag to vector to hold (relative) system ids of DiffEqs
  //!   coupled to eq among other DiffEqs of type coupledeq
  //! \tparam depvar_msg Error message key to use on missing coupled depvar
  //! \param[in] in Parser input
  //! \param[in,out] stack Grammar stack to work with
  //! \param[in] missing Error message key to use on missing coupled equation
  //!   if the coupling is required. Pass MsgKey::OPTIONAL as missing if the
  //!   coupling is optional.
  //! \details This function computes and assigns the relative system id of a
  //!   an equation coupled to another equation. The two equations being coupled
  //!   are given by the two template arguments 'eq' and 'coupledeq'. The goal
  //!   is to compute the id of the coupled eq that identifies it among
  //!   potentially multiple coupled equations. Ths id is simply an integer
  //!   which is a relative index, starting from zero, and uniquely identifies
  //!   the equation system coupled to eq. As a result, when eq is instantiated,
  //!   this id can be used to query any detail of the user configuration for
  //!   the coupledeq equation coupled to eq.
  //! \note This function is intended to be called at the very end of parsing,
  //!   when all equation systems and their configuration have been parsed so
  //!   the dependent variables, identifying equation systems, and the
  //!   specifications of their coupling, via referring to dependent variables
  //!   of equations coupled, have all been specified and thus known.
  //! \note This function is intended to be called once per coupling an equation
  //!   type (a system) to another equation (system). Example: if the velocity
  //!   eqation is coupled to three other equation systems (say, position,
  //!   dissipation, mixmassfracbeta, this function must be called three times
  //!   with eq = velocity, and coupledeq = position, dissipation, and
  //!   mixmassfracbeta so that the proper coupling information is setup for
  //!   each of the three couplings.
  template< typename eq, typename coupledeq, typename id, MsgKey depvar_msg,
            typename Input, typename Stack >
  static void couple( const Input& in, Stack& stack, MsgKey missing )
  {
    // get coupled eq configuration
    const auto& ceq = stack.template get< tag::param, eq, coupledeq >();
    // get dependent variables of coupled equation
    const auto& coupled_depvar =
      stack.template get< tag::param, coupledeq, tag::depvar >();
    // get access to coupled eq relative id vector to be filled
    auto& coupled_eq_id = stack.template get< tag::param, eq, id >();

    // Note that the size of depvar vector must equal the size of the coupledeq
    // configuration vector, ceq. The depvar stores the dependent variables
    // (characters) of all the eqs configured, e.g., potentially multiple eqs.
    // The coupledeq config vector stores the dependent variables each eq is
    // potentially coupled to. Thus the position/index in these two vectors are
    // used to identify which system is being coupled to which other system. If
    // ceq[i] = '-', then eq[i] is NOT coupled to ceq[i]. ceq must be filled
    // after a particular equation system block is finished parsing by
    // check_coupling, which does error checking but is supposed to fill in the
    // coupledeq depvar. If there is no coupling, it must put in '-'.
    Assert( ( ceq.size() ==
              stack.template get< tag::param, eq, tag::depvar >().size() ),
            "Size mismatch" );

    // Find relative system ids for all coupledeqs coupled to eqs. Note that we
    // loop through all ceqs (whose size equals to depvar, the number of eqs
    // configured) and try to find all ceq depvars among the coupledeq depvars.
    for (auto ev : ceq) {   // for all depvars coupled to eq
      std::size_t c = 0;
      for (auto cv : coupled_depvar) {  // for all depvars of a coupled eq
        if (ev == cv) coupled_eq_id.push_back( c );
        ++c;
      }
      // If eq is not coupled we put in a large number as a placeholder to
      // keep the vector sizes of ceq and coupled_eq_id equal, this id
      // will surely trigger a problem if ends up being used.
      if (ev == '-')
        coupled_eq_id.push_back( std::numeric_limits<std::size_t>::max() );
    }

    // Ensure all coupled ids are filled.
    Assert( ceq.size() == coupled_eq_id.size(), "Not all coupled eqs found" );

    // Error out if the coupling is required and at least one of the
    // potentially multiple eqs is not coupled. Since this function is called
    // for potentially multiple eqs, ceq is a vector. If the coupling between eq
    // and coupledeq is required, ceq should not contatain a '-'.
    if ( missing != MsgKey::OPTIONAL &&
         std::any_of( ceq.cbegin(), ceq.cend(),
                      []( char v ){ return v == '-'; } ) )
    {
      Message< Stack, ERROR, depvar_msg >( stack, in );
    }
  }

  //! Query if equation 'eq' has been coupled to equation 'coupledeq'
  //! \tparam eq Tag of the equation to query
  //! \tparam coupledeq Tag of the equation that is potentially coupled to
  //!   equation 'eq'
  //! \param[in] stack Grammar stack to work with
  //! \return True if equation 'eq' is coupled to equation 'coupledeq'
  //! \note Always the eq system that is parsed last is interrogated.
  template< typename eq, typename coupledeq, typename Stack >
  static bool coupled( const Stack& stack ) {
    return stack.template get< tag::param, eq, coupledeq >().back() != '-';
  }

  //! Query number of components of coupled equation
  //! \tparam eq Tag of the equation that is coupled
  //! \tparam coupledeq Tag of the equation that is coupled to equation 'eq'
  //! \tparam id Tag to access the coupled equation 'eq' (relative) ids, see
  //!   tk::grm::couple.
  //! \param[in] stack Grammar stack to work with
  //! \return Number of scalar components of coupled equation
  //! \note Always the eq system that is parsed last is interrogated.
  template< typename eq, typename coupledeq, typename id, typename Stack >
  static std::size_t ncomp_coupled( const Stack& stack ) {
    Assert( (coupled< eq, coupledeq >( stack )), "Eq must be coupled" );
    // Query relative id of coupled eq
    auto cid = stack.template get< tag::param, eq, id >().back();
    return stack.template get< tag::component, coupledeq >().at( cid );
  }

  //! Rule used to trigger action
  struct check_velocity : pegtl::success {};
  //! \brief Do error checking on the velocity eq block
  template<>
  struct action< check_velocity > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      // if there was a velocity eq block defined
      if (neq.get< tag::velocity >() > 0) {
        // Ensure a coupled position model is configured
        couple< tag::velocity,
            tag::position, tag::position_id, MsgKey::POSITION_DEPVAR >
          ( in, stack, MsgKey::OPTIONAL );
        // Compute equation id if a coupled dissipation model is configured
        couple< tag::velocity,
            tag::dissipation, tag::dissipation_id, MsgKey::DISSIPATION_DEPVAR >
          ( in, stack, MsgKey::OPTIONAL );
        // Compute equation id if a coupled mass fraction model is configured
        couple< tag::velocity, tag::mixmassfracbeta, tag::mixmassfracbeta_id,
                MsgKey::MIXMASSFRACBETA_DEPVAR >
          ( in, stack, MsgKey::OPTIONAL );
        // Error out if no dependent variable to solve for was selected
        const auto& solve =
          stack.template get< tag::param, tag::velocity, tag::solve >();
        if (solve.size() != neq.get< tag::velocity >())
          Message< Stack, ERROR, MsgKey::NOSOLVE >( stack, in );

        // Increase number of components by the number of particle densities if
        // we solve for particle momentum, coupled to mass fractions. This is
        // to allocate storage for particle velocity as variables derived from
        // momentum.
        if (!solve.empty() && solve.back() == walker::ctr::DepvarType::PRODUCT)
        {
          //! Error out if not coupled to mixmassfracbeta
          if (!coupled< tag::velocity, tag::mixmassfracbeta >( stack )) {
            Message< Stack, ERROR, MsgKey::MIXMASSFRACBETA_DEPVAR >(stack,in);
          } else {
            // access number of components of velocity eq just parsed
            auto& ncomp = stack.template get< tag::component, tag::velocity >();
            // query number of components of coupled mixmassfracbeta model
            auto nc = ncomp_coupled< tag::velocity, tag::mixmassfracbeta,
                                     tag::mixmassfracbeta_id >( stack );
            // Augment storage of velocity equation, solving for momentum by
            // the number of scalar components the coupled mixmassfracbeta mix
            // model solves for. The magic number, 4, below is
            // MixMassFractionBeta::NUMDERIVED + 1, and the 3 is the number of
            // velocity components derived from momentum.
            ncomp.back() += (nc/4)*3;
          }
        }

        // Set C0 = 2.1 if not specified
        auto& C0 = stack.template get< tag::param, tag::velocity, tag::c0 >();
        if (C0.size() != neq.get< tag::velocity >()) C0.push_back( 2.1 );
        // Set SLM if not specified
        auto& variant =
          stack.template get< tag::param, tag::velocity, tag::variant >();
        if (variant.size() != neq.get< tag::velocity >())
          variant.push_back( walker::ctr::VelocityVariantType::SLM );

        // Set gravity to {0,0,0} if unspecified
        auto& gravity =
          stack.template get< tag::param, tag::velocity, tag::gravity >();
        if (gravity.size() != neq.get< tag::velocity >())
          gravity.push_back( { 0.0, 0.0, 0.0 } );
      }
    }
  };

  //! Rule used to trigger action
  struct check_position : pegtl::success {};
  //! \brief Do error checking on the position eq block and compute coupling
  template<>
  struct action< check_position > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      // if there was a position eq block defined
      if (neq.get< tag::position >() > 0) {
        // Ensure a coupled velocity model is configured
        couple< tag::position,
            tag::velocity, tag::velocity_id, MsgKey::VELOCITY_DEPVAR >
          ( in, stack, MsgKey::VELOCITY_MISSING );
        // Error out if no dependent variable to solve for was selected
        const auto& solve =
          stack.template get< tag::param, tag::position, tag::solve >();
        if (solve.size() != neq.get< tag::position >())
          Message< Stack, ERROR, MsgKey::NOSOLVE >( stack, in );
      }
    }
  };

  //! Rule used to trigger action
  struct check_dissipation : pegtl::success {};
  //! \brief Do error checking on the dissipation eq block and compute coupling
  template<>
  struct action< check_dissipation > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      using eq = tag::dissipation;
      using param = tag::param;
      // if there was a dissipation eq block defined
      if (neq.get< eq >() > 0) {
        // Ensure a coupled velocity model is configured
        couple< eq, tag::velocity, tag::velocity_id,
                MsgKey::VELOCITY_DEPVAR >
              ( in, stack, MsgKey::VELOCITY_MISSING );
        // Set C3 if not specified
        auto& C3 = stack.template get< param, eq, tag::c3 >();
        if (C3.size() != neq.get< eq >()) C3.push_back( 1.0 );
        // Set C4 if not specified
        auto& C4 = stack.template get< param, eq, tag::c4 >();
        if (C4.size() != neq.get< eq >()) C4.push_back( 0.25 );
        // Set COM1 if not specified
        auto& COM1 = stack.template get< param, eq, tag::com1 >();
        if (COM1.size() != neq.get< eq >()) COM1.push_back( 0.44 );
        // Set COM2 if not specified
        auto& COM2 = stack.template get< param, eq, tag::com2 >();
        if (COM2.size() != neq.get< eq >()) COM2.push_back( 0.9 );
      }
    }
  };

  //! Rule used to trigger action
  struct check_mixmassfracbeta : pegtl::success {};
  //! \brief Do error checking on a mass fraction eq block
  template<>
  struct action< check_mixmassfracbeta > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using walker::deck::neq;
      using eq = tag::mixmassfracbeta;
      // Error out if no dependent variable to solve for was selected
      const auto& solve = stack.template get< tag::param, eq, tag::solve >();
      if (solve.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOSOLVE >( stack, in );
      // if there was a mixmassfracbeta eq block defined
      if (neq.get< eq >() > 0) {
        // Compute equation id if a coupled velocity model is configured
        couple< eq, tag::velocity, tag::velocity_id,
                MsgKey::VELOCITY_DEPVAR >
              ( in, stack, MsgKey::OPTIONAL );
        // Compute equation id if a coupled dissipation model is configured
        couple< eq, tag::dissipation, tag::dissipation_id,
                MsgKey::DISSIPATION_DEPVAR >
              ( in, stack, MsgKey::OPTIONAL );
      }
    }
  };

  //! Rule used to trigger action
  struct check_mixdirichlet : pegtl::success {};
  //! \brief Error checks for the mixdirichlet sde
  template<>
  struct action< check_mixdirichlet > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using eq = tag::mixdirichlet;
      using walker::deck::neq;

      // Ensure correct size for parameter vector rho
      auto& rho = stack.template get< tag::param, eq, tag::rho >().back();
      auto ncomp = stack.template get< tag::component, eq >().back();
      if (rho.size() != ncomp-2)
        Message< Stack, ERROR, MsgKey::MIXDIR_RHO >( stack, in );

      // If normalization is not set, set default
      auto& normv = stack.template get< tag::param, eq, tag::normalization >();
      if (normv.size() != neq.get< eq >())
        normv.push_back( walker::ctr::NormalizationType::LIGHT );

      // Sort parameter vector rho to correct order depending on normalization
      if (normv.back() == walker::ctr::NormalizationType::HEAVY)
        std::sort( rho.begin(), rho.end() );
      else
        std::sort( rho.begin(), rho.end(), std::greater< tk::real >() );
    }
  };

  //! Rule used to trigger action
  template< typename eq, typename coupledeq >
  struct check_coupling : pegtl::success {};
  //! Put in coupled eq depvar as '-' if no coupling is given
  template< typename eq, typename coupledeq >
  struct action< check_coupling< eq, coupledeq > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      auto& ceq = stack.template get< tag::param, eq, coupledeq >();
      const auto& depvar = stack.template get< tag::param, eq, tag::depvar >();
      // A depvar of '-' means no coupling. This keeps the coupled eq vector
      // size the same as that of the depvar vector, so we keep track of which
      // eq is coupled to which coupled eq (and also which eq is not coupled to
      // any other eq).
      if (depvar.size() > ceq.size()) ceq.push_back( '-' );
      Assert( depvar.size() == ceq.size(), "Vector size mismatch" );
    }
  };

  //! Rule used to trigger action
  struct velocity_defaults : pegtl::success {};
  //! \brief Set defaults for the velocity model
  template<>
  struct action< velocity_defaults > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      // Set number of components: always 3 velocity components
      auto& ncomp = stack.template get< tag::component, tag::velocity >();
      ncomp.push_back( 3 );
    }
  };

  //! Rule used to trigger action
  struct position_defaults : pegtl::success {};
  //! \brief Set defaults for all Lagrangian particle position models
  template<>
  struct action< position_defaults > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      // Set number of components: always 3 position components
      auto& ncomp = stack.template get< tag::component, tag::position >();
      ncomp.push_back( 3 );
      // Set RNG if no RNG has been selected (not all position models use this,
      // so don't impose on the user to define one). Pick a Random123 generator,
      // as that is always available.
      using walker::deck::neq;
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
                                               tag::init,
                                               tag::spike >,
                        // start new vector or vectors of beta parameters for a
                        // potential jointbeta initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::init,
                                               tag::betapdf >,
                        // start new vector or vectors of gamma parameters for a
                        // potential jointgamma initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::init,
                                               tag::gamma >,
                        // start new vector or vectors of gaussian parameters
                        // for a potential jointgaussian initpolicy
                        tk::grm::start_vector< tag::param,
                                               eq,
                                               tag::init,
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

  //! SDE parameter vector
  template< class keyword,
            template< class, class, class... > class check_vector,
            class eq, class param, class... xparams >
  struct sde_parameter_vector :
         tk::grm::parameter_vector< use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    check_vector,
                                    eq,
                                    param,
                                    xparams... > {};

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
                             tag::init,
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
                             tag::init,
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
                             tag::init,
                             tag::gamma > > > {};

  //! scan icdirichlet ... end block
  template< class eq >
  struct icdirichlet :
         pegtl::if_must<
           tk::grm::readkw< use< kw::icdirichlet >::pegtl_string >,
           // parse a dirichletpdf ... end block (there can be multiple)
           tk::grm::block<
             use< kw::end >,
             sde_parameter_vector< kw::dirichletpdf,
                                   tk::grm::check_dirichletpdf,
                                   eq,
                                   tag::init,
                                   tag::dirichlet > > > {};

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
                             tag::init,
                             tag::gaussian > > > {};

  //! scan icjointgaussian ... end block
  template< class eq >
  struct icjointgaussian :
         pegtl::if_must<
           tk::grm::readkw< use< kw::icjointgaussian >::pegtl_string >,
           // parse a jointgaussian ... end block
           tk::grm::block< use< kw::end >,
                           sde_parameter_vector< kw::sde_mean,
                                                 tk::grm::check_vector,
                                                 eq,
                                                 tag::init,
                                                 tag::mean >,
                           sde_parameter_vector< kw::sde_cov,
                                                 tk::grm::check_vector,
                                                 eq,
                                                 tag::init,
                                                 tag::cov > > > {};

  //! Error checks after an equation ... end block has been parsed
  template< class eq, class... extra_checks  >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           tk::grm::register_eq< eq >,
           // performe extra pegtl actions, e.g., performing extra checks
           extra_checks...,
           // do error checking on this block
           tk::grm::check_eq< eq >,
           // do error checking on the init policy
           tk::grm::check_init< eq > > {};

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
                           icjointgaussian< tag::diagou >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tk::grm::check_vector,
                                                 tag::diagou,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_theta,
                                                 tk::grm::check_vector,
                                                 tag::diagou,
                                                 tag::theta >,
                           sde_parameter_vector< kw::sde_mu,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::ou >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tk::grm::check_vector,
                                                 tag::ou,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_theta,
                                                 tk::grm::check_vector,
                                                 tag::ou,
                                                 tag::theta >,
                           sde_parameter_vector< kw::sde_mu,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::skewnormal >,
                           sde_parameter_vector< kw::sde_T,
                                                 tk::grm::check_vector,
                                                 tag::skewnormal,
                                                 tag::timescale >,
                           sde_parameter_vector< kw::sde_sigmasq,
                                                 tk::grm::check_vector,
                                                 tag::skewnormal,
                                                 tag::sigmasq >,
                           sde_parameter_vector< kw::sde_lambda,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::beta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::beta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::beta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::numfracbeta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::numfracbeta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::numfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
                                                 tag::numfracbeta,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tk::grm::check_vector,
                                                 tag::numfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_rcomma,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::massfracbeta >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::massfracbeta,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::massfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
                                                 tag::massfracbeta,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tk::grm::check_vector,
                                                 tag::massfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_r,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::mixnumfracbeta >,
                           sde_parameter_vector< kw::sde_bprime,
                                                 tk::grm::check_vector,
                                                 tag::mixnumfracbeta,
                                                 tag::bprime >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::mixnumfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappaprime,
                                                 tk::grm::check_vector,
                                                 tag::mixnumfracbeta,
                                                 tag::kappaprime >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tk::grm::check_vector,
                                                 tag::mixnumfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_rcomma,
                                                 tk::grm::check_vector,
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
                           tk::grm::policy< use,
                                            use< kw::solve >,
                                            ctr::Depvar,
                                            tag::mixmassfracbeta,
                                            tag::solve >,
                           icdelta< tag::mixmassfracbeta >,
                           icbeta< tag::mixmassfracbeta >,
                           icgamma< tag::mixmassfracbeta >,
                           icgaussian< tag::mixmassfracbeta >,
                           icjointgaussian< tag::mixmassfracbeta >,
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
                                                 tk::grm::check_vector,
                                                 tag::mixmassfracbeta,
                                                 tag::bprime >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::mixmassfracbeta,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappaprime,
                                                 tk::grm::check_vector,
                                                 tag::mixmassfracbeta,
                                                 tag::kappaprime >,
                           sde_parameter_vector< kw::sde_rho2,
                                                 tk::grm::check_vector,
                                                 tag::mixmassfracbeta,
                                                 tag::rho2 >,
                           sde_parameter_vector< kw::sde_r,
                                                 tk::grm::check_vector,
                                                 tag::mixmassfracbeta,
                                                 tag::r >,
                           sde_parameter_vector< kw::mean_gradient,
                                                 tk::grm::check_mean_gradient,
                                                 tag::mixmassfracbeta,
                                                 tag::mean_gradient >,
                           tk::grm::process<
                             use< kw::velocity >,
                             tk::grm::Store_back< tag::param,
                                                  tag::mixmassfracbeta,
                                                  tag::velocity >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::dissipation >,
                             tk::grm::Store_back< tag::param,
                                                  tag::mixmassfracbeta,
                                                  tag::dissipation >,
                             pegtl::alpha > >,
           check_errors< tag::mixmassfracbeta,
                         tk::grm::check_vector_exists<
                           tag::mixmassfracbeta,
                           tag::hydrotimescales,
                           tk::grm::MsgKey::HYDROTIMESCALES >,
                         tk::grm::check_vector_exists<
                           tag::mixmassfracbeta,
                           tag::hydroproductions,
                           tk::grm::MsgKey::HYDROPRODUCTIONS >,
                         tk::grm::check_coupling< tag::mixmassfracbeta,
                                                  tag::dissipation >,
                         tk::grm::check_coupling< tag::mixmassfracbeta,
                                                  tag::velocity > > > {};

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
                           icjointgaussian< tag::gamma >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::gamma,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::gamma,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::dirichlet >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::dirichlet,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::dirichlet,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
                                                 tag::dirichlet,
                                                 tag::kappa >
                         >,
           check_errors< tag::dirichlet > > {};

  //! MixDirichlet SDE
  struct mixdirichlet :
         pegtl::if_must<
           scan_sde< use< kw::mixdirichlet >, tag::mixdirichlet >,
           tk::grm::block< use< kw::end >,
                           tk::grm::depvar< use,
                                            tag::mixdirichlet,
                                            tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::mixdirichlet >,
                           tk::grm::rng< use,
                                         use< kw::rng >,
                                         tk::ctr::RNG,
                                         tag::mixdirichlet,
                                         tag::rng >,
                           tk::grm::policy< use,
                                            use< kw::init >,
                                            ctr::InitPolicy,
                                            tag::mixdirichlet,
                                            tag::initpolicy >,
                           tk::grm::policy< use,
                                            use< kw::coeff >,
                                            ctr::CoeffPolicy,
                                            tag::mixdirichlet,
                                            tag::coeffpolicy >,
                           tk::grm::policy< use,
                                            use< kw::normalization >,
                                            ctr::Normalization,
                                            tag::mixdirichlet,
                                            tag::normalization >,
                           icdelta< tag::mixdirichlet >,
                           icbeta< tag::mixdirichlet >,
                           icgamma< tag::mixdirichlet >,
                           icdirichlet< tag::mixdirichlet >,
                           icgaussian< tag::mixdirichlet >,
                           icjointgaussian< tag::mixdirichlet >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::mixdirichlet,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::mixdirichlet,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappaprime,
                                                 tk::grm::check_vector,
                                                 tag::mixdirichlet,
                                                 tag::kappaprime >,
                           sde_parameter_vector< kw::sde_rho,
                                                 tk::grm::check_vector,
                                                 tag::mixdirichlet,
                                                 tag::rho >
                         >,
           check_errors< tag::mixdirichlet,
                         tk::grm::check_mixdirichlet > > {};

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
                           icjointgaussian< tag::gendir >,
                           sde_parameter_vector< kw::sde_b,
                                                 tk::grm::check_vector,
                                                 tag::gendir,
                                                 tag::b >,
                           sde_parameter_vector< kw::sde_S,
                                                 tk::grm::check_vector,
                                                 tag::gendir,
                                                 tag::S >,
                           sde_parameter_vector< kw::sde_kappa,
                                                 tk::grm::check_vector,
                                                 tag::gendir,
                                                 tag::kappa >,
                           sde_parameter_vector< kw::sde_c,
                                                 tk::grm::check_vector,
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
                           icjointgaussian< tag::wrightfisher >,
                           sde_parameter_vector< kw::sde_omega,
                                                 tk::grm::check_vector,
                                                 tag::wrightfisher,
                                                 tag::omega > >,
           check_errors< tag::wrightfisher > > {};

  //! Velocity SDE
  struct velocity :
         pegtl::if_must<
           scan_sde< use< kw::velocity >, tag::velocity >,
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
                           tk::grm::policy< use,
                                            use< kw::solve >,
                                            ctr::Depvar,
                                            tag::velocity,
                                            tag::solve >,
                           tk::grm::policy< use,
                                            use< kw::variant >,
                                            ctr::VelocityVariant,
                                            tag::velocity,
                                            tag::variant >,
                           sde_parameter_vector< kw::gravity,
                                                 tk::grm::check_gravity,
                                                 tag::velocity,
                                                 tag::gravity >,
                           icdelta< tag::velocity >,
                           icbeta< tag::velocity >,
                           icgamma< tag::velocity >,
                           icgaussian< tag::velocity >,
                           icjointgaussian< tag::velocity >,
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
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::mixmassfracbeta >,
                             tk::grm::Store_back< tag::param,
                                                  tag::velocity,
                                                  tag::mixmassfracbeta >,
                             pegtl::alpha > >,
           check_errors< tag::velocity,
                         tk::grm::velocity_defaults,
                         tk::grm::check_coupling< tag::velocity,
                                                  tag::position >,
                         tk::grm::check_coupling< tag::velocity,
                                                  tag::dissipation >,
                         tk::grm::check_coupling< tag::velocity,
                                                  tag::mixmassfracbeta > > > {};

  //! position equation
  struct position :
         pegtl::if_must<
           scan_sde< use< kw::position >, tag::position >,
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
                           tk::grm::policy< use,
                                            use< kw::solve >,
                                            ctr::Depvar,
                                            tag::position,
                                            tag::solve >,
                           icdelta< tag::position >,
                           icbeta< tag::position >,
                           icgamma< tag::position >,
                           icgaussian< tag::position >,
                           icjointgaussian< tag::position >,
                           tk::grm::process<
                             use< kw::velocity >,
                             tk::grm::Store_back< tag::param,
                                                  tag::position,
                                                  tag::velocity >,
                             pegtl::alpha > >,
           check_errors< tag::position, 
                         tk::grm::position_defaults,
                         tk::grm::check_coupling< tag::position,
                                                  tag::velocity > > > {};

  //! dissipation equation
  struct dissipation :
         pegtl::if_must<
           scan_sde< use< kw::dissipation >, tag::dissipation >,
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
                           tk::grm::process<
                             use< kw::sde_c3 >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::c3 > >,
                           tk::grm::process<
                             use< kw::sde_c4 >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::c4 > >,
                           tk::grm::process<
                             use< kw::sde_com1 >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::com1 > >,
                           tk::grm::process<
                             use< kw::sde_com2 >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::com2 > >,
                           icdelta< tag::dissipation >,
                           icbeta< tag::dissipation >,
                           icgamma< tag::dissipation >,
                           icgaussian< tag::dissipation >,
                           icjointgaussian< tag::dissipation >,
                           tk::grm::process<
                             use< kw::velocity >,
                             tk::grm::Store_back< tag::param,
                                                  tag::dissipation,
                                                  tag::velocity >,
                             pegtl::alpha > >,
           check_errors< tag::dissipation,
                         tk::grm::dissipation_defaults,
                         tk::grm::check_coupling< tag::dissipation,
                                                  tag::velocity > > > {};

  //! stochastic differential equations
  struct sde :
         pegtl::sor< dirichlet,
                     mixdirichlet,
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
               tk::grm::check_dissipation,
               tk::grm::check_mixmassfracbeta >,
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
