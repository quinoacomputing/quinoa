// *****************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test suite grammar definition
  \details   Random number generator test suite input deck grammar definition.
  We use the Parsing Expression Grammar Template Library (PEGTL) to create the
  grammar and the associated parser. Word of advice: read from the bottom up.
*/
// *****************************************************************************
#ifndef RNGTestInputDeckGrammar_h
#define RNGTestInputDeckGrammar_h

#include "CommonGrammar.h"
#include "Keywords.h"
#include "QuinoaConfig.h"

#ifdef HAS_MKL
#include "MKLGrammar.h"
#endif

#ifdef HAS_RNGSSE2
  #include "RNGSSEGrammar.h"
#endif

#include "Random123Grammar.h"

namespace rngtest {

extern ctr::InputDeck g_inputdeck_defaults;

//! RNGTest input deck facilitating user input for testing RNGs
namespace deck {

  //! \brief Specialization of tk::grm::use for RNGTest's control file parser
  //! \author J. Bakosi
  template< typename keyword >
  using use = tk::grm::use< keyword,
                            ctr::InputDeck::keywords1,
                            ctr::InputDeck::keywords2,
                            ctr::InputDeck::keywords3 >;

} // ::deck
} // ::rngtest

namespace tk {
namespace grm {

  // Note that PEGTL action specializations must be in the same namespace as the
  // template being specialized. See http://stackoverflow.com/a/3052604.

  // RNGTest's InputDeck actions

  //! Rule used to trigger action
  template< class Option, typename... tags >
  struct store_rngtest_option : pegtl::success {};
  //! \brief Put option in state at position given by tags
  //! \details This is simply a wrapper around tk::grm::store_option passing the
  //!    stack defaults for rngtest.
  //! \author J. Bakosi
  template< class Option, typename... tags >
  struct action< store_rngtest_option< Option, tags... > > {
    template< typename Input, typename Stack >
    static void apply( const Input& in, Stack& stack ) {
      store_option< Stack, rngtest::deck::use, Option, rngtest::ctr::InputDeck,
                    Input, tags... >
                  ( stack, in, rngtest::g_inputdeck_defaults );
    }
  };

} // ::grm
} // ::tk

namespace rngtest {

//! RNGTest input deck facilitating user input for testing RNGs
namespace deck {

  // RNGTest's InputDeck grammar

  //! \brief Match the inside of rngs ... end block
  //! \author J. Bakosi
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

  // \brief Match TestU01 batteries
  //! \author J. Bakosi
  template< typename battery_kw >
  struct testu01 :
         pegtl::if_must<
           tk::grm::scan< typename battery_kw::pegtl_string,
                          tk::grm::store_rngtest_option< ctr::Battery,
                                                         tag::selected,
                                                         tag::battery > >,
           pegtl::sor< tk::grm::block< use< kw::end >, rngs >,
                       tk::grm::msg< tk::grm::MsgType::ERROR,
                                     tk::grm::MsgKey::UNFINISHED > > > {};

  //! \brief Match all batteries
  //! \author J. Bakosi
  struct battery :
         pegtl::sor< testu01< use< kw::smallcrush > >,
                     testu01< use< kw::crush > >,
                     testu01< use< kw::bigcrush > > > {};

  //! \brief All keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< tk::grm::title< use >, battery > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  //! \author J. Bakosi
  struct read_file :
         tk::grm::read_file< keywords, tk::grm::ignore > {};

} // deck::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
