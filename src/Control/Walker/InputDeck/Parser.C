// *****************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
// *****************************************************************************

#include <ostream>
#include <vector>
#include <type_traits>

#include "NoWarning/pegtl.h"

#include "Print.h"
#include "Tags.h"
#include "ContainerUtil.h"
#include "Walker/Types.h"
#include "Walker/InputDeck/InputDeck.h"
#include "Walker/InputDeck/Parser.h"
#include "Walker/InputDeck/Grammar.h"

namespace tk {
namespace grm {

extern tk::Print g_print;

} // grm::
} // tk::

using walker::InputDeckParser;

InputDeckParser::InputDeckParser( const tk::Print& print,
                                  const ctr::CmdLine& cmdline,
                                  ctr::InputDeck& inputdeck ) :
  FileParser( cmdline.get< tag::io, tag::control >() )
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line stack
//! \param[inout] inputdeck Input deck stack where data is stored during parsing
//! \author  J. Bakosi
// *****************************************************************************
{
  // Create InputDeck (a tagged tuple) to store parsed input
  ctr::InputDeck id( cmdline );

  // Reset parser's output stream to that of print's. This is so that mild
  // warnings emitted during parsing can be output using the pretty printer.
  // Usually, errors and warnings are simply accumulated during parsing and
  // printed during diagnostics after the parser has finished. Howver, in some
  // special cases we can provide a more user-friendly message right during
  // parsing since there is more information available to construct a more
  // sensible message. This is done in e.g., tk::grm::store_option. Resetting
  // the global g_print, to that of passed in as the constructor argument allows
  // not to have to create a new pretty printer, but use the existing one.
  tk::grm::g_print.reset( print.save() );

  // Parse input file and populate the underlying tagged tuple
  pegtl::read_parser p( m_filename );
  p.parse< deck::read_file, tk::grm::action >( id );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, id.get< tag::error >() );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  inputdeck = std::move( id );

  // Filter out repeated statistics
  tk::unique( inputdeck.get< tag::stat >() );
}
