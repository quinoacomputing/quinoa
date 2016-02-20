//******************************************************************************
/*!
  \file      src/Main/InciterPrint.C
  \author    J. Bakosi
  \date      Fri 05 Feb 2016 06:04:20 AM MST
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
//******************************************************************************

#include <boost/algorithm/string/replace.hpp>
#include <boost/iterator/iterator_traits.hpp>

#include "PDE.h"
#include "InciterPrint.h"

using inciter::InciterPrint;

namespace inciter {

extern std::vector< PDE > g_pdes;

} // inciter::

void
InciterPrint::inthead( const std::string& title,
                       const std::string& name,
                       const std::string& legend,
                       const std::string& head ) const
//******************************************************************************
//  Print time integration header
//! \param[in] title Section title
//! \param[in] name Section name
//! \param[in] legend Legend to print
//! \param[in] head Head to append
//! \author J. Bakosi
//******************************************************************************
{
  section( title, name );
  std::string l( legend );
  boost::replace_all( l, "\n", "\n" + m_item_indent );
  raw( m_item_indent + l + head );
}

void
InciterPrint::pdes( const std::string& title, const std::vector< std::vector<
  std::pair< std::string, std::string > > >& info ) const
//******************************************************************************
//  Print configuration of a stack of partial differential equations
//! \param[in] title Title to use
//! \param[in] info Info vector to use
//! \author J. Bakosi
//******************************************************************************
{
  if ( !info.empty() ) {
    std::stringstream ss;
    ss << title << " (" << g_pdes.size() << ")";
    section( ss.str() );
    for (std::size_t e=0; e<info.size(); ++e) {
      subsection( info[e][0].first );
      for (std::size_t l=1; l<info[e].size(); ++l)
        m_stream << m_item_name_value_fmt % m_item_indent
                    % info[e][l].first % info[e][l].second;
      if (e < info.size()-1) endsubsection();
    }
  }
}
