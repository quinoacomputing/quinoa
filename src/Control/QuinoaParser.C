//******************************************************************************
/*!
  \file      src/Control/QuinoaParser.C
  \author    J. Bakosi
  \date      Fri Sep 20 11:00:54 2013
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
#ifdef NDEBUG
  pegtl::dummy_parse_file< grm::read_file >( m_filename, m_control );
#else
  pegtl::basic_parse_file< grm::read_file >( m_filename, m_control );
#endif

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

void
QuinoaParser::echo() const
//******************************************************************************
//  Echo problem setup
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.section("Title", m_control.get<ctr::title>());
  m_print.item("Control file", m_control.get<ctr::io,ctr::control>());
}
