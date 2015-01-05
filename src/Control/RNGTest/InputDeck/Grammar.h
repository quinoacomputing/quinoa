//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:26:28 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Random number generator test suite grammar definition
  \details   Random number generator test suite input deck grammar definition.
  We use the Parsing Expression Grammar Template Library (PEGTL) to create the
  grammar and the associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at)
  for PEGTL. Word of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef RNGTestInputDeckGrammar_h
#define RNGTestInputDeckGrammar_h

#include <PEGTLParsed.h>
#include <Keywords.h>
#include <Grammar.h>

#ifdef HAS_MKL
#include <MKLGrammar.h>
#endif

#include <RNGSSEGrammar.h>

namespace rngtest {

extern ctr::InputDeck g_inputdeck_defaults;

//! RNGTest input deck facilitating user input for testing RNGs
namespace deck {

  //! PEGTLParsed type specialized to RNGTest's input deck parser
  //! \details PEGTLInputDeck is practically InputDeck equipped with PEGTL
  //!   location information so the location can be tracked during parsing.
  //! \author J. Bakosi
  using PEGTLInputDeck =
          tk::ctr::PEGTLParsed< ctr::InputDeck,
                                pegtl::file_input< ctr::Location >,
                                tag::cmd,
                                ctr::CmdLine >;

  //! \brief Specialization of tk::grm::use for RNGTest's control file parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3 >;

  // RNGTest's InputDeck State

  //! Everything is stored in Stack during parsing
  //! \author J. Bakosi
  using Stack = PEGTLInputDeck;

  // RNGTest's InputDeck actions

  //! \brief Pput option in state at position given by tags
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

  // RNGTest's InputDeck grammar

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

  // \brief Match TestU01 batteries
  //! \author J. Bakosi
  template< typename battery_kw >
  struct testu01 :
         pegtl::ifmust< tk::grm::scan< Stack,
                                       typename battery_kw::pegtl_string,
                                       store_option< ctr::Battery,
                                                     tag::selected,
                                                     tag::battery > >,
                        tk::grm::block< Stack, use< kw::end >, rngs > > {};

  //! \brief Match all batteries
  //! \author J. Bakosi
  struct battery :
         pegtl::sor< testu01< use< kw::smallcrush > >,
                     testu01< use< kw::crush > >,
                     testu01< use< kw::bigcrush > > > {};

  //! \brief All keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use >, battery > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  //! \author J. Bakosi
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore > {};

} // deck::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
