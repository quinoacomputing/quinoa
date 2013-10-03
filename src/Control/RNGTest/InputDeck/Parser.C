//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Thu Oct  3 17:35:34 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite input deck parser
  \details   Random number generator test suite input deck parser
*/
//******************************************************************************

#include <RNGTest/InputDeck/Parser.h>
#include <RNGTest/InputDeck/Grammar.h>

using namespace rngtest;

InputDeckParser::InputDeckParser(quinoa::Base& base) :
  FileParser(base, "blah")
//******************************************************************************
//  Constructor
//! \param[inout] base      Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  m_base.print.item("Control file", m_filename);
}

void
InputDeckParser::parse()
//******************************************************************************
//  Parse random number generator test suite control file
//! \author  J. Bakosi
//******************************************************************************
{
  // Parse: basic_parse_file() below gives debug info during parsing, use it for
  // debugging the parser itself, i.e., when modifying the grammar, otherwise,
  // use dummy_parse_file() which compiles faster
  //pegtl::dummy_parse_file< grm::read_file >( m_filename, m_base.control );
  //pegtl::basic_parse_file< grm::read_file >( m_filename, m_base.control );

  m_base.print.item("Parsed control file", "success");
}
