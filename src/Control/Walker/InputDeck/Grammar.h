// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3,
                            ctr::InputDeck::keywords4,
                            ctr::InputDeck::keywords5,
                            ctr::InputDeck::keywords6,
                            ctr::InputDeck::keywords7 >;

  // Walker's InputDeck state
 
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
} // ::deck
} // ::walker

namespace tk {
namespace grm {

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // Walker's InputDeck actions

  //! Rule used to trigger action
  template< class Option, typename... tags >
  struct store_walker_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for walker.
  //! \author J. Bakosi
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
  //! \author J. Bakosi
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
  //! \author J. Bakosi
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
  //!   required to have a specific size: ncomp/4, i.e., the number of
  //!   components in a block divided by four. This specific value is the only
  //!   way this is used at this time, thus the hard-coding of ncomp/4. However,
  //!   this could be abstracted away via a template argument if needed.
  //! \note This functor does not check existence of a vector. If the vector
  //!   does not even exist, it will throw an exception in DEBUG mode, while in
  //!   RELEASE mode it will attempt to access unallocated memory yielding a
  //!   segfault. The existence of the vector should be checked by
  //!   check_vector_exists first.
  //! \author J. Bakosi
  template< class eq, class vec >
  struct action< check_vector_size< eq, vec > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      const auto& vv = stack.template get< tag::param, eq, vec >();
      Assert( !vv.empty(), "Vector of vectors checked must not be empty" );
      const auto& v = vv.back();
      const auto& ncomp = stack.template get< tag::component, eq >().back();
      if (v.empty() || v.size() != ncomp/4)
        Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_eq : pegtl::success {};
  //! \brief Do general error checking on the differential equation block
  //! \details This is error checking that all equation types must satisfy.
  //! \author J. Bakosi
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
  //! \author J. Bakosi
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

    }
  };

} // ::grm
} // ::tk

namespace walker {

//! Walker input deck facilitating user input for integrating SDEs
namespace deck {

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
                                               tag::betapdf > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discr< use< kw::npar >, tag::npar >,
                     tk::grm::discr< use< kw::nstep >, tag::nstep >,
                     tk::grm::discr< use< kw::term >, tag::term >,
                     tk::grm::discr< use< kw::dt >, tag::dt >,
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
         pegtl::if_must<
           tk::grm::readkw< use< kw::walker >::pegtl_string >,
           pegtl::sor<
             tk::grm::block<
               use< kw::end >,
               discretization_parameters,
               sde,
               tk::grm::rngblock< use, rngs >,
               tk::grm::statistics< use, tk::grm::store_walker_option >,
               tk::grm::pdfs< use, tk::grm::store_walker_option > >,
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
