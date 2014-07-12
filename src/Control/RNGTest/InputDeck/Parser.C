//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:52:48 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************

#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/InputDeck/Grammar.h>

using rngtest::InputDeckParser;

InputDeckParser::InputDeckParser( const tk::Print& print,
                                  ctr::CmdLine& cmdline,
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

  // Strip input deck (and its underlying tagged tuple) from PEGTL instruments
  // and transfer it out
  inputdeck = id;

  // If we got here, parser succeeded
  print.item("Parsed control file", "success");
}
