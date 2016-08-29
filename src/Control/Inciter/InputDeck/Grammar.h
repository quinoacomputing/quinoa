// *****************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Mon 29 Aug 2016 02:10:56 PM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter's input deck grammar definition
  \details   Inciter's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef InciterInputDeckGrammar_h
#define InciterInputDeckGrammar_h

#include "CommonGrammar.h"
#include "PEGTLParsed.h"
#include "Keywords.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/InputDeck/InputDeck.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck_defaults;

//! Inciter input deck facilitating user input for computing shock hydrodynamics
namespace deck {

  //! \brief PEGTLParsed type specialized to Inciter's input deck parser
  //! \details PEGTLInputDeck is practically InputDeck equipped with PEGTL
  //!   location information so the location can be tracked during parsing.
  //! \author J. Bakosi
  using PEGTLInputDeck =
    tk::ctr::PEGTLParsed< ctr::InputDeck,
                          pegtl::file_input< ctr::Location >,
                          tag::cmd,
                          ctr::CmdLine >;

  //! \brief Specialization of tk::grm::use for Inciter's input deck parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3 >;

  // Inciter's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  //! \author J. Bakosi
  static tk::tuple::tagged_tuple< tag::transport, std::size_t,
                                  tag::poisson,   std::size_t,
                                  tag::compflow,  std::size_t > neq;

  // Inciter's InputDeck actions

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
  //! \details This is error checking that all equation types must satisfy. For
  //!   more specific equations, such as compressible flow, a more specialized
  //!   equation checker does and can do better error checking.
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
      // Error out if no test problem has been selected
      const auto& problem = stack.get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOINIT >
                        ( stack, value );
    }
  };

  //! \brief Set defaults and do error checking on the transport equation block
  //! \details This is error checking that only the transport equation block
  //!   must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a transport...end block has
  //!   just finished.
  //! \author J. Bakosi
  template< class eq >
  struct check_transport : pegtl::action_base< check_transport< eq > > {
    static void apply( const std::string& value, Stack& stack ) {
      // If no number of components has been selected, default to 1
      const auto& ncomp = stack.get< tag::component, eq >();
      if (ncomp.empty() || ncomp.size() != neq.get< eq >())
        stack.get< tag::component, eq >().push_back( 1 );
      // If physics type is not given, default to 'advection'
      auto& physics = stack.get< tag::param, eq, tag::advection >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( ctr::PhysicsType::ADVECTION );
      // If problem type is not given, error out
      auto& problem = stack.get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOPROBLEM >
                        ( stack, value );
    }
  };

  //! \brief Set defaults and do error checking on the compressible flow
  //!   equation block
  //! \details This is error checking that only the compressible flow equation
  //!   block must satisfy. Besides error checking we also set defaults here as
  //!   this block is called when parsing of a compflow...end block has
  //!   just finished.
  //! \author J. Bakosi
  template< class eq >
  struct check_compflow : pegtl::action_base< check_compflow< eq > > {
    static void apply( const std::string& value, Stack& stack ) {
      // Set number of components to 5 (mass, 3 x mom, energy)
      stack.get< tag::component, eq >().push_back( 5 );
      // If physics type is not given, default to 'euler'
      auto& physics = stack.get< tag::param, eq, tag::physics >();
      if (physics.empty() || physics.size() != neq.get< eq >())
        physics.push_back( ctr::PhysicsType::EULER );
      // If problem type is not given, default to 'user_defined'
      auto& problem = stack.get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        problem.push_back( ctr::ProblemType::USER_DEFINED );
      else if (problem.back() == ctr::ProblemType::VORTICAL_FLOW) {
        const auto& alpha = stack.get< tag::param, eq, tag::alpha >();
        const auto& beta = stack.get< tag::param, eq, tag::beta >();
        const auto& p0 = stack.get< tag::param, eq, tag::p0 >();
        if ( alpha.size() != problem.size() ||
             beta.size() != problem.size() ||
             p0.size() != problem.size() )
          tk::grm::Message< Stack, tk::grm::ERROR,
                            tk::grm::MsgKey::VORTICAL_UNFINISHED >
                          ( stack, value );
      }
    }
  };

  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults.
  //! \author J. Bakosi
  template< class Option, typename... tags >
  struct store_option : pegtl::action_base< store_option< Option, tags... > > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::grm::store_option< Stack, use, Option, ctr::InputDeck, tags... >
                           ( stack, value, g_inputdeck_defaults );
    }
  };

  // Inciter's InputDeck grammar

  //! scan and store_back equation keyword and option
  template< typename keyword, class eq >
  struct scan_eq :
         tk::grm::scan< Stack,
                        typename keyword::pegtl_string,
                        tk::grm::store_back_option< Stack,
                                                    use,
                                                    ctr::PDE,
                                                    tag::selected,
                                                    tag::pde > > {};

  //! Error checks after an equation...end block has been parsed
  template< class eq, template< class > class eqchecker >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           pegtl::apply< register_eq< eq > >,
           // do error checking on this block
           pegtl::apply< eqchecker< eq > > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discr< Stack, use< kw::nstep >, tag::nstep >,
                     tk::grm::discr< Stack, use< kw::term >, tag::term >,
                     tk::grm::discr< Stack, use< kw::t0 >, tag::t0 >,
                     tk::grm::discr< Stack, use< kw::dt >, tag::dt >,
                     tk::grm::interval< Stack, use< kw::ttyi >, tag::tty > > {};

  //! PDE parameter vector
  template< class keyword, class eq, class param >
  struct pde_parameter_vector :
         tk::grm::parameter_vector< Stack,
                                    use,
                                    use< keyword >,
                                    tk::grm::Store_back_back,
                                    tk::grm::start_vector,
                                    tk::grm::check_vector,
                                    eq,
                                    param > {};

  //! Dirichlet boundary conditions block
  template< class eq, class param >
  struct bc_dirichlet :
           pegtl::ifmust<
             tk::grm::readkw< Stack, use< kw::bc_dirichlet >::pegtl_string >,
             tk::grm::block<
               Stack,
               use< kw::end >,
               tk::grm::parameter_vector< Stack,
                                          use,
                                          use< kw::sideset >,
                                          tk::grm::Store_back_back,
                                          tk::grm::start_vector,
                                          tk::grm::check_vector,
                                          eq,
                                          param > > > {};

  //! initial conditions block for compressible flow
  template< class eq, class param >
  struct ic_compflow :
           pegtl::ifmust<
             tk::grm::readkw< Stack, use< kw::ic >::pegtl_string >,
             tk::grm::block<
               Stack,
               use< kw::end >,
               tk::grm::parameter_vector< Stack,
                                          use,
                                          use< kw::velocity >,
                                          tk::grm::Store_back_back,
                                          tk::grm::start_vector,
                                          tk::grm::check_vector,
                                          eq,
                                          param > > > {};

  //! put in material property for equation matching keyword
  template< typename eq, typename keyword, typename property >
  struct material_property :
         tk::grm::process< Stack, use< keyword >,
           tk::grm::Store_back< Stack, tag::param, eq, property > > {};

  //! Material properties block for compressible flow
  template< class eq >
  struct material_properties :
           pegtl::ifmust<
             tk::grm::readkw< Stack, use< kw::material >::pegtl_string >,
             tk::grm::block< Stack,
                             use< kw::end >,
                             material_property< eq, kw::id, tag::id >,
                             material_property< eq, kw::mat_gamma, tag::gamma >,
                             material_property< eq, kw::mat_mu, tag::mu >,
                             material_property< eq, kw::mat_cv, tag::cv >,
                             material_property< eq, kw::mat_k, tag::k > > > {};

  //! put in PDE parameter for equation matching keyword
  template< typename eq, typename keyword, typename p >
  struct parameter :
         tk::grm::process< Stack, use< keyword >,
           tk::grm::Store_back< Stack, tag::param, eq, p > > {};

  //! transport equation for scalars
  struct transport :
         pegtl::ifmust<
           scan_eq< use< kw::transport >, tag::transport >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::transport,
                                            tag::physics >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::transport,
                                            tag::problem >,
                          tk::grm::depvar< Stack,
                                           use,
                                           tag::transport,
                                           tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::transport >,
                           pde_parameter_vector< kw::pde_diffusivity,
                                                 tag::transport,
                                                 tag::diffusivity >,
                           pde_parameter_vector< kw::pde_lambda,
                                                 tag::transport,
                                                 tag::lambda >,
                           pde_parameter_vector< kw::pde_u0,
                                                 tag::transport,
                                                 tag::u0 > >,
           check_errors< tag::transport, check_eq > > {};

  //! Poisson partial differential equation for a scalar
  struct poisson :
         pegtl::ifmust<
           scan_eq< use< kw::poisson >, tag::poisson >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::poisson,
                                            tag::problem >,
                          tk::grm::depvar< Stack,
                                           use,
                                           tag::poisson,
                                           tag::depvar >,
                           tk::grm::component< Stack,
                                               use< kw::ncomp >,
                                               tag::poisson >,
                           bc_dirichlet< tag::poisson, tag::bc_dirichlet > >,
           check_errors< tag::poisson, check_eq > > {};

  //! compressible flow
  struct compflow :
         pegtl::ifmust<
           scan_eq< use< kw::compflow >, tag::compflow >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::physics >,
                                            ctr::Physics,
                                            tag::compflow,
                                            tag::physics >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::compflow,
                                            tag::problem >,
                           //ic_compflow< tag::compflow, tag::ic > >,
                           material_properties< tag::compflow >,
                           parameter< tag::compflow, kw::npar, tag::npar >,
                           parameter< tag::compflow, kw::pde_alpha, tag::alpha >,
                           parameter< tag::compflow, kw::pde_beta, tag::beta >,
                           parameter< tag::compflow, kw::pde_p0, tag::p0 >,
                           bc_dirichlet< tag::compflow, tag::bc_dirichlet > >,
           check_errors< tag::compflow, check_compflow > > {};

  //! partitioning ... end block
  struct partitioning :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::partitioning >::pegtl_string >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::process<
                             Stack,
                             use< kw::algorithm >,
                             store_option< tk::ctr::PartitioningAlgorithm,
                                           tag::selected,
                                           tag::partitioner >,
                             pegtl::alpha > > > {};

  //! equation types
  struct equations :
         pegtl::sor< transport, poisson, compflow > {};

  //! plotvar ... end block
  struct plotvar :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::plotvar >::pegtl_string >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::interval< Stack,
                                              use< kw::interval >,
                                              tag::field > > > {};

  //! 'inciter' block
  struct inciter :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::inciter >::pegtl_string >,
           pegtl::sor< tk::grm::block< Stack,
                         use< kw::end >,
                         discretization_parameters,
                         equations,
                         partitioning,
                         plotvar,
                         tk::grm::diagnostics< Stack, use, store_option > >,
                       pegtl::apply<
                          tk::grm::error< Stack,
                                          tk::grm::MsgKey::UNFINISHED > > > > {};

  //! \brief All keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use >, inciter > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  //! \author J. Bakosi
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore< Stack > > {};

} // deck::
} // inciter::

#endif // InciterInputDeckGrammar_h
