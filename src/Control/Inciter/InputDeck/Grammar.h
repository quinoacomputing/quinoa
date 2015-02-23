//******************************************************************************
/*!
  \file      src/Control/Inciter/InputDeck/Grammar.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:57:50 AM MST
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

  // Inciter's InputDeck actions

  // Inciter's InputDeck grammar

  //! \brief All keywords
  //! \author J. Bakosi
  struct keywords :
         pegtl::sor< tk::grm::title< Stack, use > > {};

  //! \brief Grammar entry point: parse keywords and ignores until eof
  //! \author J. Bakosi
  struct read_file :
         tk::grm::read_file< Stack, keywords, tk::grm::ignore > {};

} // deck::
} // inciter::

#endif // InciterInputDeckGrammar_h
