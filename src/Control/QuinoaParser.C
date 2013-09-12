//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Thu Sep 12 10:49:16 2013
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
  // Parse
  pegtl::dummy_parse_file<grammar::read_file>(m_filename, m_control);

  // Filter out repeated statistics
  using namespace control;
  unique(const_cast<std::vector<Product>&>(m_control.get<stats>()));
}

void
QuinoaParser::unique(std::vector<control::Product>& statistics)
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

void
QuinoaParser::echo() const
//******************************************************************************
//  Echo parsed information from quinoa control
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.section("Problem title", m_control.get<control::title>());
}
