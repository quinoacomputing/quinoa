//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Sun 06 Oct 2013 03:30:19 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's comamnd line parser
  \details   RNGTest's comamnd line parser
*/
//******************************************************************************

#include <RNGTest/CmdLine/Parser.h>
#include <RNGTest/CmdLine/Grammar.h>

using namespace rngtest;

void
CmdLineParser::parse()
//******************************************************************************
//  Parse command line file
//! \author  J. Bakosi
//******************************************************************************
{
  // Parse: basic_parse_string() below gives debug info during parsing, use it
  // for debugging the parser itself, i.e., when modifying the grammar,
  // otherwise, use dummy_parse_string() which compiles faster
  pegtl::dummy_parse_string< cmd::read_string >( m_string, m_base.control );
  //pegtl::basic_parse_string< cmd::read_string >( m_string, m_base.control );

  m_base.print.item("Parsed command line", "success");

  // Make sure mandatory arguments are set
  ErrChk(!(m_base.control.get<ctr::io, ctr::control>().empty()),
         quinoa::ExceptType::FATAL,
         "Mandatory control file not specified. "
         "Use '--control <filename>' or '-c <filename>'.");
}
