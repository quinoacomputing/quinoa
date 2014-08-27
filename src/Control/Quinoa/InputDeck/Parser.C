//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Tue 26 Aug 2014 06:26:44 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************

#include <make_unique.h>

#include <Quinoa/InputDeck/Parser.h>
#include <Quinoa/InputDeck/Grammar.h>

using quinoa::InputDeckParser;

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

  // Parse input file by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< deck::read_file >( input, id );

  // Echo errors and warnings accumulated during parsing
  echoErrors( print, id.get< tk::tag::error >() );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  inputdeck = id;

  // Filter out repeated statistics
  unique( inputdeck.get< tag::stat >() );

  // If we got here, parser succeeded
  print.item( "Parsed control file", "success" );
}

void
InputDeckParser::unique( std::vector< ctr::Product >& statistics )
//******************************************************************************
//  Make requested statistics unique
//! \param[in,out]  statistics  Vector of statistics
//! \author  J. Bakosi
//******************************************************************************
{
  std::sort( begin(statistics), end(statistics) );
  auto it = std::unique( begin(statistics), end(statistics) );
  statistics.resize( std::distance( begin(statistics), it ) );
}
