//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Parser.C
  \author    J. Bakosi
  \date      Thu Oct  3 10:19:22 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck file parser
  \details   Quinoa's input deck file parser
*/
//******************************************************************************

#include <algorithm>

#include <pegtl.hh>

#include <Quinoa/InputDeck/Parser.h>
#include <Quinoa/InputDeck/Grammar.h>

using namespace quinoa;

InputDeckParser::InputDeckParser(const std::string& filename, Base& base) :
  FileParser(filename, base)
//******************************************************************************
//  Constructor
//! \param[in]    filename  Control file name to read from
//! \param[inout] base      Essentials
//! \author  J. Bakosi
//******************************************************************************
{
  m_base.control.set<ctr::io,ctr::control>(filename);
  m_base.print.item("Control file", filename);
}

void
InputDeckParser::parse()
//******************************************************************************
//  Parse quinoa control file
//! \author  J. Bakosi
//******************************************************************************
{
  // Parse: basic_parse_file() below gives debug info during parsing, use it for
  // debugging the parser itself, i.e., when modifying the grammar, otherwise,
  // use dummy_parse_file() which compiles faster
  pegtl::dummy_parse_file< grm::read_file >( m_filename, m_base.control );
  //pegtl::basic_parse_file< grm::read_file >( m_filename, m_base.control );

  m_base.print.item("Parsed control file", "success");

  // Filter out repeated statistics
  unique(m_base.control.get<ctr::stat>());
}

void
InputDeckParser::unique(std::vector<ctr::Product>& statistics)
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
