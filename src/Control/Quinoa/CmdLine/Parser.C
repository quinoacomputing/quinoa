//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/Parser.C
  \author    J. Bakosi
  \date      Wed Oct  2 16:40:15 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's comamnd line parser
  \details   Quinoa's comamnd line parser
*/
//******************************************************************************

#include <pegtl.hh>

#include <Quinoa/CmdLine/Parser.h>
#include <Quinoa/CmdLine/Grammar.h>

using namespace quinoa;

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
  //pegtl::basic_parse_string< grm::read_string >( m_string, m_base.control );

  m_base.print.item("Parsed command line", "success");
}
