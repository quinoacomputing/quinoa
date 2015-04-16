//******************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Sun 15 Mar 2015 06:43:55 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter's input deck grammar definition
  \details   Inciter's input deck grammar definition. We use the Parsing
  Expression Grammar Template Library (PEGTL) to create the grammar and the
  associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at) for PEGTL. Word
  of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef InciterInputDeckGrammar_h
#define InciterInputDeckGrammar_h

#include <Grammar.h>
#include <PEGTLParsed.h>
#include <Keywords.h>

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
                            ctr::InputDeck::keywords >;

  // Inciter's InputDeck state

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  //! \brief Number of registered equations
  //! \details Counts the number of parsed equation blocks during parsing.
  //! \author J. Bakosi
  tk::tuple::tagged_tuple< tag::scalar, std::size_t > neq;

  // Inciter's InputDeck actions

  //! \brief Register differential equation after parsing its block
  //! \details This is used by the error checking functors (check_*) during
  //!    parsing to identify the recently-parsed block.
  //! \author J. Bakosi
  template< class eq >
  struct register_eq : pegtl::action_base< register_eq< eq > > {
    static void apply( const std::string& value, Stack& stack ) {
      ++neq.get< eq >();
    }
  };

  //! \brief Do general error checking on the differential equation block
  //! \details This is error checking that all equation types must satisfy.
  //! \author J. Bakosi
  template< class eq >
  struct check_eq : pegtl::action_base< check_eq< eq > > {
    static void apply( const std::string& value, Stack& stack ) {
      // Error out if no test problem has been selected
      const auto& problem = stack.get< tag::param, eq, tag::problem >();
      if (problem.empty() || problem.size() != neq.get< eq >())
        tk::grm::Message< Stack, tk::grm::ERROR, tk::grm::MsgKey::NOINIT >
                        ( stack, value );
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
                                                    ctr::DiffEq,
                                                    tag::selected,
                                                    tag::diffeq > > {};

  //! Error checks after an equation...end block has been parsed
  template< class eq >
  struct check_errors :
         pegtl::seq<
           // register differential equation block
           pegtl::apply< register_eq< eq > >,
           // do error checking on this block
           pegtl::apply< check_eq< eq > > > {};

  //! Discretization parameters
  struct discretization_parameters :
         pegtl::sor< tk::grm::discr< Stack, use< kw::nstep >, tag::nstep >,
                     tk::grm::discr< Stack, use< kw::term >, tag::term >,
                     tk::grm::discr< Stack, use< kw::dt >, tag::dt >,
                     tk::grm::interval< Stack, use< kw::ttyi >, tag::tty > > {};

  //! scalar transport equation
  struct scalar :
         pegtl::ifmust<
           scan_eq< use< kw::scalar >, tag::scalar >,
           tk::grm::block< Stack,
                           use< kw::end >,
                           tk::grm::policy< Stack,
                                            use,
                                            use< kw::problem >,
                                            ctr::Problem,
                                            tag::scalar,
                                            tag::problem > >,
           check_errors< tag::scalar > > {};

  //! equation types
  struct equations :
         pegtl::sor< scalar > {};

  //! 'inciter' block
  struct inciter :
         pegtl::ifmust<
           tk::grm::readkw< Stack, use< kw::inciter >::pegtl_string >,
           pegtl::sor< tk::grm::block< Stack,
                                       use< kw::end >,
                                       discretization_parameters,
                                       equations >,
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
