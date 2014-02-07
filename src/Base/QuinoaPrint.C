//******************************************************************************
/*!
  \file      src/Base/QuinoaPrint.C
  \author    J. Bakosi
  \date      Fri 07 Feb 2014 09:35:58 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     QuinoaPrint
  \details   QuinoaPrint
*/
//******************************************************************************

#include <QuinoaPrint.h>
#include <Quinoa/Options/SDE.h>

using quinoa::QuinoaPrint;

void
QuinoaPrint::RequestedStats(const std::string& msg) const
//******************************************************************************
//  Echo requested statistics if differs from default
//! \author J. Bakosi
//******************************************************************************
{
  if (m_ctr.get<tag::stat>() != ctr::InputDeckDefaults.get<tag::stat>()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : m_ctr.get<tag::stat>()) {
      m_stream <<= v;
    }
    m_stream << '\n';
  }
}

void
QuinoaPrint::EstimatedStats(const std::string& msg) const
//******************************************************************************
//  Echo estimated statistics if differs from default
//! \author J. Bakosi
//******************************************************************************
{
  if (m_ctr.get<tag::stat>() != ctr::InputDeckDefaults.get<tag::stat>()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : m_ctr.get<tag::stat>()) {
       m_stream << v;
    }
    m_stream << '\n';
  }
}

std::vector< std::string >
QuinoaPrint::SDEPolicyNames( const ctr::SDEKey& key ) const
//******************************************************************************
//  Return SDE policies names
//! \author J. Bakosi
//******************************************************************************
{
  std::vector< std::string > names;
  names.push_back( ctr::InitPolicy().name( key.get<tag::initpolicy>() ) );
  names.push_back( ctr::CoeffPolicy().name( key.get<tag::coeffpolicy>() ) );
  return names;
}
