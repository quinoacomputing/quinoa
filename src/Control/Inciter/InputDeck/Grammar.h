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
#include "ContainerUtil.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  //! \brief Specialization of tk::grm::use for Inciter's input deck parser
  template< typename keyword >
  using use = tk::grm::use< keyword, ctr::InputDeck::keywords >;

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
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );
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
      // Error out if no dependent variable has been selected
      auto& depvar = stack.template get< tag::param, eq, tag::depvar >();
      if (depvar.empty() || depvar.size() != neq.get< eq >())
        Message< Stack, ERROR, MsgKey::NODEPVAR >( stack, in );
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
  struct check_inciter : pegtl::success {};
  //! \brief Do error checking on the inciter block
  template<> struct action< check_inciter > {
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
      // if MatCG is configured, turn on reordering
      if (stack.template get< tag::discr, tag::scheme >() ==
           inciter::ctr::SchemeType::MatCG)
        stack.template get< tag::discr, tag::reorder >() = true;
      // if DGP1 is configured, set ndofs to be 4
      if (stack.template get< tag::discr, tag::scheme >() ==
           inciter::ctr::SchemeType::DGP1)
        stack.template get< tag::discr, tag::ndof >() = 4;
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

  //! Rule used to trigger action
  struct compute_refvar_idx : pegtl::success {};
  //! Compute indices of refinement variables
  //! \details This functor computes the indices in the unknown vector for all
  //!   refinement variables in the system of systems of dependent variables
  //!   after the refvar...end block has been parsed in the amr...end block.
  //!   After basic error checking, the vector at stack.get<tag::amr,tag::id>()
  //!   is filled.
  template<>
  struct action< compute_refvar_idx > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // reference variables just parsed by refvar...end block
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      // get ncomponents object from this input deck
      const auto& ncomps = stack.template get< tag::component >();
      // compute offset map associating offsets to dependent variables
      auto offsetmap = ncomps.offsetmap( stack );
      // compute number of components associated to dependent variabels
      auto ncompmap = ncomps.ncompmap( stack );
      // reference variable index vector to fill
      auto& refidx = stack.template get< tag::amr, tag::id >();
      // Compute indices for all refvars
      for (const auto& v : refvar) {    // for all reference variables parsed
        // depvar is the first char of a refvar
        auto depvar = v[0];
        // the field ID is optional and is the rest of the depvar string
        std::size_t f = (v.size()>1 ? std::stoul(v.substr(1)) : 1) - 1;
        // field ID must be less than or equal to the number of scalar
        // components configured for the eq system for this dependent variable
        if (f >= tk::cref_find( ncompmap, depvar ))
          Message< Stack, ERROR, MsgKey::NOSUCHCOMPONENT >( stack, in );
        // get offset for depvar
        auto eqsys_offset = tk::cref_find( offsetmap, depvar );
        // the index is the eq offset + field ID
        auto idx = eqsys_offset + f;
        // save refvar index in system of all systems
        refidx.push_back( idx );
      }
    }
  };

  //! Rule used to trigger action
  struct check_amr_errors : pegtl::success {};
  //! Do error checking for the amr...end block
  //! \details This is error checking that only the amr...end block
  //!   must satisfy. Besides error checking this can also set defaults
  //!   as this block is called when parsing of a amr...end block has
  //!   just finished.
  template<>
  struct action< check_amr_errors > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      // Error out if refvar size does not equal refidx size (programmer error)
      Assert( (stack.template get< tag::amr, tag::refvar >().size() ==
               stack.template get< tag::amr, tag::id >().size()),
              "The size of refvar and refidx vectors must equal" );
      const auto& initref = stack.template get< tag::amr, tag::init >();
      const auto& edgeref = stack.template get< tag::amr, tag::edge >();
      const auto& refvar = stack.template get< tag::amr, tag::refvar >();
      // Error out if initref edge list is not divisible by 2 (user error)
      if (edgeref.size() % 2 == 1)
        Message< Stack, ERROR, MsgKey::T0REFODD >( stack, in );
      // Error out if initial AMR will be a no-op (user error)
      if ( stack.template get< tag::amr, tag::t0ref >() &&
           initref.empty() && edgeref.empty() )
        Message< Stack, ERROR, MsgKey::T0REFNOOP >( stack, in );
      // Error out if timestepping AMR will be a no-op (user error)
      if ( stack.template get< tag::amr, tag::dtref >() && refvar.empty() )
        Message< Stack, ERROR, MsgKey::DTREFNOOP >( stack, in );
      // Error out if mesh refinement frequency is zero (programmer error)
      Assert( (stack.template get< tag::amr, tag::dtfreq >() > 0),
              "Mesh refinement frequency must be positive" );
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

  //! Match discretization option
  template< template< class > class use, class keyword, class Option,
            class Tag >
  struct discroption :
         tk::grm::process< use< keyword >,
                           tk::grm::store_inciter_option<
                             Option, tag::discr, Tag >,
                           pegtl::alpha > {};

  //! Discretization parameters
  struct discretization :
         pegtl::sor<
           tk::grm::discrparam< use, kw::nstep, tag::nstep >,
           tk::grm::discrparam< use, kw::term, tag::term >,
           tk::grm::discrparam< use, kw::t0, tag::t0 >,
           tk::grm::discrparam< use, kw::dt, tag::dt >,
           tk::grm::discrparam< use, kw::cfl, tag::cfl >,
           tk::grm::discrparam< use, kw::ctau, tag::ctau >,
           tk::grm::process< use< kw::fct >, 
                             tk::grm::Store< tag::discr, tag::fct >,
                             pegtl::alpha >,
           tk::grm::process< use< kw::reorder >,
                             tk::grm::Store< tag::discr, tag::reorder >,
                             pegtl::alpha >,
           tk::grm::interval< use< kw::ttyi >, tag::tty >,
           discroption< use, kw::scheme, inciter::ctr::Scheme, tag::scheme >,
           discroption< use, kw::flux, inciter::ctr::Flux, tag::flux >,
           discroption< use, kw::limiter, inciter::ctr::Limiter, tag::limiter >,
           tk::grm::discrparam< use, kw::cweight, tag::cweight >
         > {};

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

  //! initref ... end block
  struct initref :
         tk::grm::vector< use< kw::amr_initref >,
                          tk::grm::Store_back< tag::amr, tag::edge >,
                          use< kw::end >,
                          tk::grm::check_vector< tag::amr, tag::edge > > {};

  //! xminus configuring coordinate-based edge tagging for mesh refinement
  template< typename keyword, typename Tag >
  struct half_world :
         tk::grm::control< use< keyword >, pegtl::digit, tag::amr, Tag > {};

  //! coordref ... end block
  struct coordref :
           pegtl::if_must<
             tk::grm::readkw< use< kw::amr_coordref >::pegtl_string >,
             tk::grm::block< use< kw::end >,
                         half_world< kw::amr_xminus, tag::xminus >,
                         half_world< kw::amr_xplus, tag::xplus >,
                         half_world< kw::amr_yminus, tag::yminus >,
                         half_world< kw::amr_yplus, tag::yplus >,
                         half_world< kw::amr_zminus, tag::zminus >,
                         half_world< kw::amr_zplus, tag::zplus > > > {};

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
                           bc< kw::bc_outlet, tag::transport, tag::bcoutlet >,
                           bc< kw::bc_extrapolate, tag::transport,
                               tag::bcextrapolate > >,
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
                           tk::grm::depvar< use,
                                            tag::compflow,
                                            tag::depvar >,
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
                           bc< kw::bc_outlet, tag::compflow, tag::bcoutlet >,
                           bc< kw::bc_extrapolate, tag::compflow,
                               tag::bcextrapolate > >,
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

  //! equation types
  struct equations :
         pegtl::sor< transport, compflow > {};

  //! refinement variable(s) (refvar) ... end block
  struct refvars :
         pegtl::if_must<
           tk::grm::vector< use< kw::amr_refvar >,
                            tk::grm::match_depvar<
                              tk::grm::Store_back< tag::amr, tag::refvar > >,
                            use< kw::end >,
                            tk::grm::check_vector< tag::amr, tag::refvar >,
                            tk::grm::fieldvar< pegtl::alpha > >,
           tk::grm::compute_refvar_idx > {};

  //! adaptive mesh refinement (AMR) amr...end block
  struct amr :
         pegtl::if_must<
           tk::grm::readkw< use< kw::amr >::pegtl_string >,
           tk::grm::enable_amr, // enable AMR if amr...end block encountered
           tk::grm::block< use< kw::end >,
                           refvars,
                           initref,
                           coordref,
                           tk::grm::process<
                             use< kw::amr_initial >,
                             tk::grm::store_back_option< use,
                                                         ctr::AMRInitial,
                                                         tag::amr,
                                                         tag::init >,
                             pegtl::alpha >,
                           tk::grm::process<
                             use< kw::amr_error >,
                             tk::grm::store_inciter_option<
                               ctr::AMRError,
                               tag::amr, tag::error >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_t0ref >,
                             tk::grm::Store< tag::amr, tag::t0ref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtref >,
                             tk::grm::Store< tag::amr, tag::dtref >,
                             pegtl::alpha >,
                           tk::grm::process< use< kw::amr_dtfreq >,
                             tk::grm::Store< tag::amr, tag::dtfreq >,
                             pegtl::digit >
                         >,
           tk::grm::check_amr_errors > {};

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
                           discretization,
                           equations,
                           amr,
                           partitioning,
                           plotvar,
                           tk::grm::diagnostics<
                             use,
                             tk::grm::store_inciter_option > >,
                         tk::grm::check_inciter >,
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
