//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Mon 09 Sep 2013 09:53:22 PM MDT
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
  using namespace control;
  std::cout << "Parsed from " << m_filename << ":\n" << std::setfill('-')
            << std::setw(13+m_filename.length()) << "-" << std::endl;

  std::cout << " * Title: " << m_control.get<title>() << std::endl;
}
