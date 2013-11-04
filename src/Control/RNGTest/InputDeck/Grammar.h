//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Sun 03 Nov 2013 02:31:57 PM MST
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
#include <PEGTLParsed.h>
#include <RNGTest/Types.h>
#include <RNGTest/InputDeck/Keywords.h>

namespace rngtest {
namespace deck {

  using namespace pegtl;
  using namespace tk::grm;

  //! PEGTLParsed type specialized to RNGTest's input deck parser
  using PEGTLInputDeck = quinoa::ctr::PEGTLParsed< ctr::InputDeck,
                                                   file_input< ctr::Location >,
                                                   ctr::cmd,
                                                   ctr::CmdLine >;

  // RNGTest's InputDeck State

  //! Everything is stored in Stack during parsing
  using Stack = PEGTLInputDeck;

  // RNGTest's InputDeck actions

  //! put option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_option : action_base< store_option<OptionType, tags...> > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::Option<OptionType> opt;
      //! Emit warning on overwrite
      if (stack.get< tags... >() != ctr::InputDeckDefaults.get< tags... >()) {
        std::cout << "\n>>> PARSER WARNING: Multiple definitions for '"
                  << opt.group() << "' option. Overwriting '"
                  << opt.name(stack.get<tags...>()) << "' with '"
                  << opt.name(opt.value(value)) << "'.\n\n";
      }
      stack.set< tags... >( opt.value( value ) );
    }
  };

  //! push back option in state at position given by tags
  template< class OptionType, typename... tags >
  struct store_back_option : action_base< store_option<OptionType, tags...> > {
    static void apply( const std::string& value, Stack& stack ) {
      tk::Option< OptionType > opt;
      stack.push_back< tags... >( opt.value( value ) );
    }
  };

  // RNGTest's InputDeck grammar

  //! title
  struct title :
         ifmust< readkw< kw::title::pegtl_string >,
                 quoted< Stack, Set< Stack, ctr::title > > > {};

  //! rngs block
  struct rngs :
         process_rng< Stack,
                      store_back_option< quinoa::ctr::RNG,
                                         ctr::selected, ctr::rng >,
                      kw::seed::pegtl_string,
                      Store_back< Stack, ctr::param, ctr::rng, ctr::seed > > {};

  //! mkl_uniform_method
  struct mkl_uniform_method :
         process< Stack,
                  kw::mkl_uniform_method::pegtl_string,
                  store_option< ctr::MKLUniformMethod,
                                ctr::selected, ctr::mkl_uniform_method >,
                  alpha > {};

  // smallcrush block
  struct smallcrush :
         ifmust< scan< kw::smallcrush::pegtl_string,
                       store_option< ctr::Battery,
                                     ctr::selected,
                                     ctr::battery > >,
                 block< Stack,
                        rngs,
                        mkl_uniform_method > > {};

  // crush block
  struct crush :
         ifmust< scan< kw::crush::pegtl_string,
                       store_option< ctr::Battery,
                                     ctr::selected,
                                     ctr::battery > >,
                 block< Stack, rngs > > {};

  // bigcrush block
  struct bigcrush :
         ifmust< scan< kw::bigcrush::pegtl_string,
                       store_option< ctr::Battery,
                                     ctr::selected,
                                     ctr::battery > >,
                 block< Stack, rngs > > {};

  //! batteries
  struct battery :
         sor< smallcrush,
              crush,
              bigcrush > {};

  //! main keywords
  struct keywords :
         sor< title,
              battery > {};

  //! ignore: comments and empty lines
  struct ignore :
         sor< comment, until<eol, space> > {};

  //! entry point: parse keywords and ignores until eof
  struct read_file :
         until< eof, sor<keywords, ignore, unknown<Stack,Error::KEYWORD>> > {};

} // deck::
} // rngtest::

#endif // RNGTestInputDeckGrammar_h
