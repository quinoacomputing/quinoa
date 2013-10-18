//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Fri Oct 18 12:24:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************

#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/InputDeck/Grammar.h>

using namespace rngtest;

InputDeckParser::InputDeckParser(const tk::Print& print,
                                 const std::unique_ptr< ctr::CmdLine >& cmdline,
                                 std::unique_ptr< ctr::InputDeck >& inputdeck) :
  FileParser( cmdline->get<ctr::io, ctr::control>() )
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
  std::unique_ptr< deck::PEGTLInputDeck > pid( new deck::PEGTLInputDeck(input) );

  // Parse input file by populating the underlying tagged tuple:
  // basic_parse() below gives debug info during parsing, use it for debugging
  // the parser itself, i.e., when modifying the grammar, otherwise, use
  // dummy_parse() to compile faster
  pegtl::dummy_parse< deck::read_file >( input, *pid );

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // by creating a unique_ptr to the base class (InputDeck) and transfer it out
  inputdeck = std::unique_ptr< ctr::InputDeck >( std::move(pid) );

  // If we got here, parser succeeded
  print.item("Parsed control file", "success");
}
