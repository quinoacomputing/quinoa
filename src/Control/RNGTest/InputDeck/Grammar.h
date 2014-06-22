//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Wed 11 Jun 2014 01:37:40 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite grammar definition
  \details   Random number generator test suite input deck grammar definition.
  We use the Parsing Expression Grammar Template Library (PEGTL) to create the
  grammar and the associated parser. Credit goes to Colin Hirsch (pegtl@cohi.at)
  for PEGTL. Word of advice: read from the bottom up.
*/
//******************************************************************************
#ifndef RNGTestInputDeckGrammar_h
#define RNGTestInputDeckGrammar_h

#include <Macro.h>
#include <Exception.h>
#include <Option.h>
#include <PEGTLParsed.h>
#include <RNGTest/Types.h>
#include <RNGTest/InputDeck/Keywords.h>
#include <Grammar.h>

#ifdef HAS_MKL
#include <MKLGrammar.h>
#endif

#include <RNGSSEGrammar.h>

namespace rngtest {

extern ctr::InputDeck g_inputdeck_defaults;

namespace deck {

  //! PEGTLParsed type specialized to RNGTest's input deck parser
  using PEGTLInputDeck =
          quinoa::ctr::PEGTLParsed< ctr::InputDeck,
                                    pegtl::file_input< ctr::Location >,
                                    tag::cmd,
                                    ctr::CmdLine >;

  // RNGTest's InputDeck State

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  // RNGTest's InputDeck actions

  //! put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : pegtl::action_base< store_option<OptionType,tags...> > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::grm::store_option< Stack, OptionType, ctr::InputDeck, tags... >
                           ( stack, value, g_inputdeck_defaults );
    }
  };

  // RNGTest's InputDeck grammar

  //! title
  struct title :
         pegtl::ifmust< tk::grm::readkw< tk::kw::title::pegtl_string >,
                        tk::grm::quoted< Stack,
                                         tk::grm::Set<Stack, tag::title> > > {};

  //! rngs
  struct rngs :
         pegtl::sor<
                     #ifdef HAS_MKL
                     tk::mkl::rngs< Stack,
                                    tag::selected, tk::tag::rng,
                                    tag::param, tk::tag::rngmkl >,
                     #endif
                     tk::rngsse::rngs< Stack,
                                       tag::selected, tk::tag::rng,
                                       tag::param, tk::tag::rngsse > > {};

  // TestU01 batteries
  template< typename battery_kw >
  struct testu01 :
         pegtl::ifmust< tk::grm::scan< typename battery_kw::pegtl_string,
                                       store_option< ctr::Battery,
                                                     tag::selected,
                                                     tag::battery > >,
                        tk::grm::block< Stack, rngs > > {};

  //! batteries
  struct battery :
         pegtl::sor< testu01< kw::smallcrush >,
                     testu01< kw::crush >,
                     testu01< kw::bigcrush > > {};

  //! main keywords
  struct keywords :
         pegtl::sor< title,
                     battery > {};

  //! ignore: comments and empty lines
  struct ignore :
         pegtl::sor< tk::grm::comment,
                     pegtl::until< pegtl::eol, pegtl::space > > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         tk::grm::read_file< Stack, keywords, ignore > {};

} // deck::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
