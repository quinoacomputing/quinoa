// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Inciter's input deck grammar definition
  \details   Inciter's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterInputDeckGrammar_h
#define InciterInputDeckGrammar_h

#include <limits>
#include <cmath>

#include "CommonGrammar.h"
#include "Keywords.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  //! \brief Specialization of tk::grm::use for Inciter's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3,
                            ctr::InputDeck::keywords4,
                            ctr::InputDeck::keywords5 >;

  // Inciter's InputDeck state

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  static tk::tuple::tagged_tuple< tag::transport, std::size_t,
                                  tag::compflow,  std::size_t > neq;

} // ::deck
} // ::inciter

namespace tk {
namespace grm {

  using namespace tao;

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // Inciter's InputDeck actions

  //! Rule used to trigger action
  template< class eq > struct register_inciter_eq : pegtl::success {};
  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  template< class eq >
  struct action< register_inciter_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& ) {
      using inciter::deck::neq;
      ++neq.get< eq >();
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_inciter_eq : pegtl::success {};
  //! \brief Do general error checking on the differential equation block
  //! \details This is error checking that all equation types must satisfy. For
  //!   more specific equations, such as compressible flow, a more specialized
  //!   equation checker does and can do better error checking.
  template< class eq >
  struct action< check_inciter_eq< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      // Error out if no dependent variable has been selected
      const auto& depvar = stack.template get< tag::param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );
      // Error out if no number of components has been selected
      const auto& ncomp = stack.template get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NONCOMP >( stack, in );
      // Error out if no test problem has been selected
      const auto& problem = stack.template get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOINIT >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_transport : pegtl::success {};
  //! \brief Set defaults and do error checking on the transport equation block
  //! \details This is error checking that only the transport equation block
  //!   must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a transport...end block has
  //!   just finished.
  template< class eq >
  struct action< check_transport< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      // Error out if no dependent variable has been selected
      auto& depvar = stack.template get< tag::param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        depvar.push_back( 'c' );
      // If no number of components has been selected, default to 1
      auto& ncomp = stack.template get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        ncomp.push_back( 1 );
      // If physics type is not given, default to 'advection'
      auto& physics = stack.template get< tag::param, eq, tag::physics >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( inciter::ctr::PhysicsType::ADVECTION );
      // If physics type is advection-diffusion, check for correct number of
      // advection velocity, shear, and diffusion coefficients
      if (physics.back() == inciter::ctr::PhysicsType::ADVDIFF) {
        auto& u0 = stack.template get< tag::param, eq, tag::u0 >();
        if (u0.back().size() != ncomp.back())  // must define 1 component
          Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
        auto& diff = stack.template get< tag::param, eq, tag::diffusivity >();
        if (diff.back().size() != 3*ncomp.back())  // must define 3 components
          Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
        auto& lambda = stack.template get< tag::param, eq, tag::lambda >();
        if (lambda.back().size() != 2*ncomp.back()) // must define 2 shear comps
          Message< Stack, ERROR, MsgKey::WRONGSIZE >( stack, in );
      }
      // If problem type is not given, error out
      auto& problem = stack.template get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NOPROBLEM >( stack, in );
      // Error check Dirichlet boundary condition block for all transport eq
      // configurations
      for (const auto& s : stack.template get< tag::param, eq, tag::bcdir >())
        if (s.empty())
          Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class eq > struct check_compflow : pegtl::success {};
  //! \brief Set defaults and do error checking on the compressible flow
  //!   equation block
  //! \details This is error checking that only the compressible flow equation
  //!   block must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a compflow...end block has
  //!   just finished.
  template< class eq >
  struct action< check_compflow< eq > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      // Set number of components to 5 (mass, 3 x mom, energy)
      stack.template get< tag::component, eq >().push_back( 5 );
      // If physics type is not given, default to 'euler'
      auto& physics = stack.template get< tag::param, eq, tag::physics >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( inciter::ctr::PhysicsType::EULER );
      // If problem type is not given, default to 'user_defined'
      auto& problem = stack.template get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        problem.push_back( inciter::ctr::ProblemType::USER_DEFINED );
      else if (problem.back() == inciter::ctr::ProblemType::VORTICAL_FLOW) {
        const auto& alpha = stack.template get< tag::param, eq, tag::alpha >();
        const auto& beta = stack.template get< tag::param, eq, tag::beta >();
        const auto& p0 = stack.template get< tag::param, eq, tag::p0 >();
        if ( alpha.size() != problem.size() ||
             beta.size() != problem.size() ||
             p0.size() != problem.size() )
          Message< Stack, ERROR, MsgKey::VORTICAL_UNFINISHED >( stack, in );
      }
      else if (problem.back() == inciter::ctr::ProblemType::NL_ENERGY_GROWTH) {
        const auto& alpha = stack.template get< tag::param, eq, tag::alpha >();
        const auto& betax = stack.template get< tag::param, eq, tag::betax >();
        const auto& betay = stack.template get< tag::param, eq, tag::betay >();
        const auto& betaz = stack.template get< tag::param, eq, tag::betaz >();
        const auto& kappa = stack.template get< tag::param, eq, tag::kappa >();
        const auto& r0 = stack.template get< tag::param, eq, tag::r0 >();
        const auto& ce = stack.template get< tag::param, eq, tag::ce >();
        if ( alpha.size() != problem.size() ||
             betax.size() != problem.size() ||
             betay.size() != problem.size() ||
             betaz.size() != problem.size() ||
             kappa.size() != problem.size() ||
             r0.size() != problem.size() ||
             ce.size() != problem.size() )
          Message< Stack, ERROR, MsgKey::ENERGY_UNFINISHED >( stack, in);
      }
      else if (problem.back() == inciter::ctr::ProblemType::RAYLEIGH_TAYLOR) {
        const auto& alpha = stack.template get< tag::param, eq, tag::alpha >();
        const auto& betax = stack.template get< tag::param, eq, tag::betax >();
        const auto& betay = stack.template get< tag::param, eq, tag::betay >();
        const auto& betaz = stack.template get< tag::param, eq, tag::betaz >();
        const auto& kappa = stack.template get< tag::param, eq, tag::kappa >();
        const auto& p0 = stack.template get< tag::param, eq, tag::p0 >();
        const auto& r0 = stack.template get< tag::param, eq, tag::r0 >();
        if ( alpha.size() != problem.size() ||
             betax.size() != problem.size() ||
             betay.size() != problem.size() ||
             betaz.size() != problem.size() ||
             kappa.size() != problem.size() ||
             p0.size() != problem.size() ||
             r0.size() != problem.size() )
          Message< Stack, ERROR, MsgKey::RT_UNFINISHED >( stack, in);
      }
      // Error check Dirichlet boundary condition block for all compflow
      // configurations
      for (const auto& s : stack.template get< tag::param, eq, tag::bcdir >())
        if (s.empty())
          Message< Stack, ERROR, MsgKey::BC_EMPTY >( stack, in );
    }
  };

  //! Rule used to trigger action
  template< class Option, typename...tags >
  struct store_inciter_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for inciter.
  template< class Option, typename... tags >
  struct action< store_inciter_option< Option, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      store_option< Stack, inciter::deck::use, Option, inciter::ctr::InputDeck,
                    Input, tags... >
                  ( stack, in, inciter::g_inputdeck_defaults );
    }
  };

  //! Rule used to trigger action
  struct check_dt : pegtl::success {};
  //! \brief Do error checking on setting the time step calculation policy
  template<> struct action< check_dt > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      using inciter::deck::neq;
      using inciter::g_inputdeck_defaults;
      // Error out if no dt policy has been selected
      const auto& dt = stack.template get< tag::discr, tag::dt >();
      const auto& cfl = stack.template get< tag::discr, tag::cfl >();
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) <
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) <
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, ERROR, MsgKey::NODT >( stack, in );
      // If both dt and cfl are given, warn that dt wins over cfl
      if ( std::abs(dt - g_inputdeck_defaults.get< tag::discr, tag::dt >()) >
            std::numeric_limits< tk::real >::epsilon() &&
          std::abs(cfl - g_inputdeck_defaults.get< tag::discr, tag::cfl >()) >
            std::numeric_limits< tk::real >::epsilon() )
        Message< Stack, WARNING, MsgKey::MULDT >( stack, in );
    }
  };

  //! Rule used to trigger action
  struct enable_amr : pegtl::success {};
  //! Enable adaptive mesh refinement (AMR)
  template<>
  struct action< enable_amr > {
    template< typename Input, typename Stack >
    static void apply( const Input&, Stack& stack ) {
      stack.template get< tag::amr, tag::amr >() = true;
    }
  };

} // ::grm
} // ::tk

namespace inciter {

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  using namespace tao;

  // Inciter's InputDeck grammar

  //! scan and store_back equation keyword and option
  template< typename keyword, class eq >
  struct scan_eq :
         tk::grm::scan< typename keyword::pegtl_string,
                        tk::grm::store_back_option< use,
                                                    ctr::PDE,
                                                    tag::selected,
                                                    tag::pde > > {};

  //! Error checks after an equation...end block has been parsed
  template< class eq, template< class > class eqchecker >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           tk::grm::register_inciter_eq< eq >,
           // do error checking on this block
           eqchecker< eq > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discr< use< kw::nstep >, tag::nstep >,
                     tk::grm::discr< use< kw::term >, tag::term >,
                     tk::grm::discr< use< kw::t0 >, tag::t0 >,
                     tk::grm::discr< use< kw::dt >, tag::dt >,
                     tk::grm::discr< use< kw::cfl >, tag::cfl >,
                     tk::grm::discr< use< kw::ctau >, tag::ctau >,
                     tk::grm::interval< use< kw::ttyi >, tag::tty > > {};

  //! PDE parameter vector
  template< class keyword, class eq, class param >
  struct pde_parameter_vector :
         tk::grm::parameter_vector< use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq,
                                    param > {};

  //! Boundary conditions block
  template< class keyword, class eq, class param >
  struct bc :
         pegtl::if_must<
           tk::grm::readkw< typename use< keyword >::pegtl_string >,
           tk::grm::block<
             use< kw::end >,
             tk::grm::parameter_vector< use,
                                        use< kw::sideset >,
                                        tk::grm::Store_back_back,
                                        tk::grm::start_vector,
                                        tk::grm::check_vector,
                                        eq,
                                        param > > > {};

  //! initial conditions block for compressible flow
  template< class eq, class param >
  struct ic_compflow :
           pegtl::if_must<
             tk::grm::readkw< use< kw::ic >::pegtl_string >,
             tk::grm::block<
               use< kw::end >,
               tk::grm::parameter_vector< use,
                                          use< kw::velocity >,
                                          tk::grm::Store_back_back,
                                          tk::grm::start_vector,
                                          tk::grm::check_vector,
                                          eq,
                                          param > > > {};

  //! put in material property for equation matching keyword
  template< typename eq, typename keyword, typename property >
  struct material_property :
         tk::grm::process< use< keyword >,
           tk::grm::Store_back< tag::param, eq, property > > {};

  //! Material properties block for compressible flow
  template< class eq >
  struct material_properties :
           pegtl::if_must<
             tk::grm::readkw< use< kw::material >::pegtl_string >,
             tk::grm::block< use< kw::end >,
                             material_property< eq, kw::id, tag::id >,
                             material_property< eq, kw::mat_gamma, tag::gamma >,
                             material_property< eq, kw::mat_mu, tag::mu >,
                             material_property< eq, kw::mat_cv, tag::cv >,
                             material_property< eq, kw::mat_k, tag::k > > > {};

  //! put in PDE parameter for equation matching keyword
  template< typename eq, typename keyword, typename p,
            class kw_type = tk::grm::number >
  struct parameter :
         tk::grm::process< use< keyword >,
                           tk::grm::Store_back< tag::param, eq, p >,
                           kw_type > {};

  //! transport equation for scalars
  struct transport :
         pegtl::if_must<
           scan_eq< use< kw::transport >, tag::transport >,
           tk::grm::block< use< kw::end >,
                           tk::grm::policy< use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::transport,
                                            tag::physics >,
                           tk::grm::policy< use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::transport,
                                            tag::problem >,
                          tk::grm::depvar< use,
                                           tag::transport,
                                           tag::depvar >,
                           tk::grm::component< use< kw::ncomp >,
                                               tag::transport >,
                           pde_parameter_vector< kw::pde_diffusivity,
                                                 tag::transport,
                                                 tag::diffusivity >,
                           pde_parameter_vector< kw::pde_lambda,
                                                 tag::transport,
                                                 tag::lambda >,
                           pde_parameter_vector< kw::pde_u0,
                                                 tag::transport,
                                                 tag::u0 >,
                           bc< kw::bc_dirichlet, tag::transport, tag::bcdir >,
                           bc< kw::bc_sym, tag::transport, tag::bcsym >,
                           bc< kw::bc_inlet, tag::transport, tag::bcinlet >,
                           bc< kw::bc_outlet, tag::transport, tag::bcoutlet > >,
           check_errors< tag::transport, tk::grm::check_transport > > {};

  //! compressible flow
  struct compflow :
         pegtl::if_must<
           scan_eq< use< kw::compflow >, tag::compflow >,
           tk::grm::block< use< kw::end >,
                           tk::grm::policy< use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::compflow,
                                            tag::physics >,
                           tk::grm::policy< use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::compflow,
                                            tag::problem >,
                           //ic_compflow< tag::compflow, tag::ic > >,
                           material_properties< tag::compflow >,
                           parameter< tag::compflow, kw::npar, tag::npar,
                                      pegtl::digit >,
                           parameter< tag::compflow, kw::pde_alpha, tag::alpha >,
                           parameter< tag::compflow, kw::pde_p0, tag::p0 >,
                           parameter< tag::compflow, kw::pde_betax, tag::betax >,
                           parameter< tag::compflow, kw::pde_betay, tag::betay >,
                           parameter< tag::compflow, kw::pde_betaz, tag::betaz >,
                           parameter< tag::compflow, kw::pde_beta, tag::beta >,
                           parameter< tag::compflow, kw::pde_r0, tag::r0 >,
                           parameter< tag::compflow, kw::pde_ce, tag::ce >,
                           parameter< tag::compflow, kw::pde_kappa, tag::kappa >,
                           bc< kw::bc_dirichlet, tag::compflow, tag::bcdir >,
                           bc< kw::bc_sym, tag::compflow, tag::bcsym >,
                           bc< kw::bc_inlet, tag::compflow, tag::bcinlet >,
                           bc< kw::bc_outlet, tag::compflow, tag::bcoutlet > >,
           check_errors< tag::compflow, tk::grm::check_compflow > > {};

  //! partitioning ... end block
  struct partitioning :
         pegtl::if_must<
           tk::grm::readkw< use< kw::partitioning >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::process<
                             use< kw::algorithm >,
                             tk::grm::store_inciter_option<
                               tk::ctr::PartitioningAlgorithm,
                               tag::selected,
                               tag::partitioner >,
                             pegtl::alpha > > > {};

  //! discretization ... end block
  struct discretization :
         pegtl::if_must<
           tk::grm::readkw< use< kw::discretization >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::process<
                             use< kw::scheme >,
                             tk::grm::store_inciter_option<
                               inciter::ctr::Scheme,
                               tag::selected,
                               tag::scheme >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::fct >,
                             tk::grm::Store< tag::discr, tag::fct >,
                             pegtl::alpha > > > {};

  //! equation types
  struct equations :
         pegtl::sor< transport, compflow > {};

  //! adaptive mesh refinement (AMR) amr...end block
  struct amr :
         pegtl::if_must<
           tk::grm::readkw< use< kw::amr >::pegtl_string >,
           tk::grm::enable_amr, // enable AMR if amr...end block encountered
           tk::grm::block< use< kw::end >,
                           tk::grm::process<
                             use< kw::amr_initial >,
                             tk::grm::store_back_option< use,
                                                         ctr::AMRInitial,
                                                         tag::amr,
                                                         tag::init >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::amr_uniform_levels >,
                             tk::grm::Store< tag::amr, tag::levels >,
                             pegtl::digit >,
                           tk::grm::process<
                             use< kw::amr_error >,
                             tk::grm::store_inciter_option<
                               ctr::AMRError,
                               tag::amr, tag::error >,
                             pegtl::alpha > > > {};

  //! plotvar ... end block
  struct plotvar :
         pegtl::if_must<
           tk::grm::readkw< use< kw::plotvar >::pegtl_string >,
           tk::grm::block< use< kw::end >,
                           tk::grm::process< use< kw::filetype >,
                                             tk::grm::store_inciter_option<
                                               tk::ctr::FieldFile,
                                               tag::selected,
                                               tag::filetype >,
                                             pegtl::alpha >,
                           tk::grm::interval< use< kw::interval >,
                                              tag::field > > > {};

  //! 'inciter' block
  struct inciter :
         pegtl::if_must<
           tk::grm::readkw< use< kw::inciter >::pegtl_string >,
           pegtl::sor<
             pegtl::seq< tk::grm::block<
                           use< kw::end >,
                           discretization_parameters,
                           equations,
                           amr,
                           partitioning,
                           discretization,
                           plotvar,
                           tk::grm::diagnostics<
                             use,
                             tk::grm::store_inciter_option > >,
                         tk::grm::check_dt >,
            tk::grm::msg< tk::grm::MsgType::ERROR,
                          tk::grm::MsgKey::UNFINISHED > > > {};

  //! \brief All keywords
  struct keywords :
         pegtl::sor< tk::grm::title< use >, inciter > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< keywords, tk::grm::ignore > {};

} // deck::
} // inciter::

#endif // InciterInputDeckGrammar_h
