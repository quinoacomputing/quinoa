//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 11:03:16 PM MDT
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
#include <Grammar.h>
#include <RNGTest/InputDeck/Types.h>
#include <RNGTest/InputDeck/Keywords.h>

namespace rngtest {
//! RNGTest's grammar definition: state, actions, grammar
namespace grm {

  using namespace pegtl;
  using namespace quinoa;
  using namespace quinoa::grm;

  // State

  //! Everything is stored in Stack during parsing
  using Stack = InputDeck;

  // Actions

  //! put value in state at position given by tags without conversion
  template< typename... tags >
  struct set : action_base< set<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.set<tags...>(value);
    }
  };

  //! put value in state at position given by tags
  template< typename... tags >
  struct store : action_base< store<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store<tags...>(value);
    }
  };

  //! convert and push back value to vector in state at position given by tags
  template< typename...tags >
  struct store_back : action_base< store_back<tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      stack.store_back<tags...>(value);
    }
  };

  //! convert and put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : action_base< store_option<OptionType, tags...> > {
    static void apply(const std::string& value, Stack& stack) {
      quinoa::ctr::Option<OptionType> opt;
      //! Emit warning on overwrite
      if (stack.get<tags...>() != RNGTestDefaults.get<tags...>()) {
        std::cout << "\n>>> PARSER WARNING: Multiple definitions for '"
                  << opt.group() << "' option. Overwriting '"
                  << opt.name(stack.get<tags...>()) << "' with '"
                  << opt.name(opt.value(value)) << "'.\n\n";
      }
      stack.set<tags...>(opt.value(value));
    }
  };

  // Grammar

  //! rng: one of the random number generators
  struct rng :
         sor< quinoa::kw::mkl_mcg31::pegtl_string,
              quinoa::kw::mkl_r250::pegtl_string,
              quinoa::kw::mkl_mrg32k3a::pegtl_string,
              quinoa::kw::mkl_mcg59::pegtl_string,
              quinoa::kw::mkl_wh::pegtl_string,
              quinoa::kw::mkl_mt19937::pegtl_string,
              quinoa::kw::mkl_mt2203::pegtl_string,
              quinoa::kw::mkl_sfmt19937::pegtl_string,
              quinoa::kw::mkl_sobol::pegtl_string,
              quinoa::kw::mkl_niederr::pegtl_string,
              quinoa::kw::mkl_iabstract::pegtl_string,
              quinoa::kw::mkl_dabstract::pegtl_string,
              quinoa::kw::mkl_sabstract::pegtl_string,
              quinoa::kw::mkl_nondeterm::pegtl_string > {};

  // common to all RNG test suites
  struct rngtest_common :
         quinoa::grm::list< Stack,
                            kw::end::pegtl_string,
                            kw::rngs::pegtl_string,
                            store_back<ctr::generator> > {};

//   //! title
//   struct title :
//          ifmust< read<kw::title::pegtl_string>,
//                  quoted<Stack,set<ctr::title>> > {};

  // smallcrush block
  struct smallcrush :
         ifmust< quinoa::grm::parse<kw::smallcrush::pegtl_string,
                                    store_option<sel::Battery,
                                                 ctr::selected,
                                                 ctr::battery>>,
                 rngtest_common > {};

  // crush block
  struct crush :
         ifmust< quinoa::grm::parse<kw::crush::pegtl_string,
                                    store_option<sel::Battery,
                                                 ctr::selected,
                                                 ctr::battery>>,
                 rngtest_common > {};

  // bigcrush block
  struct bigcrush :
         ifmust< quinoa::grm::parse<kw::bigcrush::pegtl_string,
                                    store_option<sel::Battery,
                                                 ctr::selected,
                                                 ctr::battery>>,
                 rngtest_common > {};

  //! batteries
  struct battery :
         sor< smallcrush,
              crush,
              bigcrush > {};

  //! main keywords
  struct keywords :
         sor< success /*title,
              battery*/ > {};

  //! ignore: comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore, unknown<Stack,Error::KEYWORD>> > {};

} // grm::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
