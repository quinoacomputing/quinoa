//******************************************************************************
/*!
  \file      src/Control/Walker/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Fri 16 Jan 2015 06:05:07 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Walker's input deck file parser
  \details   Walker's input deck file parser
*/
//******************************************************************************

#include <Walker/InputDeck/Parser.h>
#include <Walker/InputDeck/Grammar.h>

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
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  print.item("Control file", m_filename);

  // Create PEGTL file input from std::string
  pegtl::file_input< ctr::Location > input( m_filename );

  // Create PEGTLInputDeck to store parsed input deck data which derives from
  // InputDeck and has location() used during parse
  deck::PEGTLInputDeck id( input, cmdline );

  // Reset parser's output stream to that of print's
  tk::grm::g_print.reset( print.save() );

  // Parse input file by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< deck::read_file >( input, id );

  // Echo errors and warnings accumulated during parsing
  diagnostics( print, id.get< tag::error >() );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  inputdeck = std::move( id );

  // Filter out repeated statistics
  tk::ctr::unique( inputdeck.get< tag::stat >() );

  // If we got here, the parser has succeeded
  print.item( "Parsed control file", "success" );
}
