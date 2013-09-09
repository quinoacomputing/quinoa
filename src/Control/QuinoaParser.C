//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Mon Sep  9 08:58:26 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa control file parser
  \details   Quinoa control file parser
*/
//******************************************************************************

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
  using namespace control;

  //std::cout << "==== PARSE START ====" << std::endl;
#ifdef NDEBUG
  //pegtl::dummy_parse_file<grammar::read_file>(m_filename, m_control);
#else  // NDEBUG
  //pegtl::basic_parse_file<grammar::read_file>(m_filename, m_control);
#endif // NDEBUG
  //std::cout << "==== PARSE END ====" << std::endl << std::endl;

  // Filter out repeated statistics
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
