//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Sat 25 Jan 2014 02:31:37 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************

#include <make_unique.h>

#include <Quinoa/InputDeck/Parser.h>
#include <Quinoa/InputDeck/Grammar.h>

using quinoa::InputDeckParser;

InputDeckParser::InputDeckParser(const tk::Print& print,
                                 std::unique_ptr< ctr::CmdLine > cmdline,
                                 std::unique_ptr< ctr::InputDeck >& inputdeck) :
  FileParser( cmdline->get<tag::io, tag::control>() )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  print.item("Control file", m_filename);

  // Create PEGTL file input from std::string
  pegtl::file_input< ctr::Location > input( m_filename );

  // Create std::unique_ptr behind which to store parsed input deck data:
  // PEGTLInputDeck derives from InputDeck and has location() used during parse
  std::unique_ptr< deck::PEGTLInputDeck >
    pid( tk::make_unique< deck::PEGTLInputDeck >( input, *cmdline ) );

  // Parse input file by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< deck::read_file >( input, *pid );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // by creating a unique_ptr to the base class (InputDeck) and transfer it out
  inputdeck = std::unique_ptr< ctr::InputDeck >( std::move(pid) );

  // Filter out repeated statistics
  unique( inputdeck->get< tag::stat >() );

  // If we got here, parser succeeded
  print.item("Parsed control file", "success");
}

void
InputDeckParser::unique( std::vector< ctr::Product >& statistics )
//******************************************************************************
//  Make requested statistics unique
//! \param[in,out]  statistics  Vector of statistics
//! \author  J. Bakosi
//******************************************************************************
{
  std::sort(statistics.begin(), statistics.end());
  auto it = std::unique(statistics.begin(), statistics.end());
  statistics.resize(std::distance(statistics.begin(), it));
}
