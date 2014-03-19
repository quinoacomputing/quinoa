//******************************************************************************
/*!
  \file      src/Base/Print.C
  \author    J. Bakosi
  \date      Wed Mar 19 10:30:51 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Print
  \details   Print
*/
//******************************************************************************

#include <Print.h>

using tk::Print;

void
Print::header(const std::string& title) const
//******************************************************************************
//  Print header
//! \author J. Bakosi
//******************************************************************************
{
  m_stream << m_header_fmt % boost::io::group(std::setfill('='), "");
  m_stream << std::endl;
  m_stream << m_header_fmt % title;
  m_stream << std::endl;
  m_stream << m_header_fmt % boost::io::group(std::setfill('='), "");
}

void
Print::part(const std::string& title) const
//******************************************************************************
//  Print part header: title
//! \author J. Bakosi
//******************************************************************************
{
  using std::operator+;
  std::string::size_type half_length = title.size()/2 + 1;
  std::string s(half_length, '-');
  std::string underline(s + " o " + s);
  std::string upper(title);
  std::transform(title.begin(), title.end(), upper.begin(), ::toupper);
  upper = "< " + upper + " >";
  m_stream << m_part_fmt % upper;
  m_stream << m_part_underline_fmt % underline;
}
