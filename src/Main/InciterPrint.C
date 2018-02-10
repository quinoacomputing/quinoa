// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.C
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************

#include "NoWarning/replace.h"
#include <boost/iterator/iterator_traits.hpp>

#include "CGPDE.h"
#include "InciterPrint.h"

using inciter::InciterPrint;

namespace inciter {

extern std::vector< CGPDE > g_cgpde;

} // inciter::

void
InciterPrint::inthead( const std::string& t,
                       const std::string& name,
                       const std::string& legend,
                       const std::string& head ) const
// *****************************************************************************
//  Print time integration header
//! \param[in] t Section title
//! \param[in] name Section name
//! \param[in] legend Legend to print
//! \param[in] head Head to append
// *****************************************************************************
{
  section( t, name );
  std::string l( legend );
  boost::replace_all( l, "\n", "\n" + m_item_indent );
  raw( m_item_indent + l + head );
}

void
InciterPrint::pdes( const std::string& t, const std::vector< std::vector<
  std::pair< std::string, std::string > > >& info ) const
// *****************************************************************************
//  Print configuration of a stack of partial differential equations
//! \param[in] t Title to use
//! \param[in] info Info vector to use
// *****************************************************************************
{
  if ( !info.empty() ) {
    std::stringstream ss;
    ss << t << " (" << g_cgpde.size() << ")";
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
