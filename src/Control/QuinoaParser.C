//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Thu 26 Sep 2013 08:36:59 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa control file parser
  \details   Quinoa control file parser
*/
//******************************************************************************

#include <algorithm>

#include <pegtl.hh>

#include <QuinoaParser.h>
#include <QuinoaGrammar.h>
#include <Control.h>

using namespace quinoa;

void
QuinoaParser::parse()
//******************************************************************************
//  Parse quinoa control file
//! \author  J. Bakosi
//******************************************************************************
{
  // Parse: basic_parse_file() below gives debug info during parsing, use it for
  // debugging the parser itself, i.e., when modifying the grammar, otherwise,
  // use dummy_parse_file() which compiles faster
  pegtl::dummy_parse_file< grm::read_file >( m_filename, m_control );
  //pegtl::basic_parse_file< grm::read_file >( m_filename, m_control );

  m_print.item("Parsed control file", "success");
  m_print.endpart();

  // Filter out repeated statistics
  unique(const_cast<std::vector<ctr::Product>&>(m_control.get<ctr::stats>()));
}

void
QuinoaParser::unique(std::vector<ctr::Product>& statistics)
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
