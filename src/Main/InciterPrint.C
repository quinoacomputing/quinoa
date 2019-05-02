// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.C
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************

#include <regex>

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
  l = std::regex_replace( l, std::regex("\n"), "\n" + m_item_indent );
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

void InciterPrint::refvar( const std::vector< std::string >& rvar,
                           const std::vector< std::size_t >& refidx )
// *****************************************************************************
// Print mesh refinement variables and their indices in the unknown vector
//! \param[in] rvar Refinement variable name list
//! \param[in] refidx Refinement variable index (location in data array) list
// *****************************************************************************
{
  Assert( rvar.size() == refidx.size(), "Size mismatch" );

  if (rvar.empty()) return;

  std::string c;
  for (std::size_t i=0; i<rvar.size(); ++i)
    c += rvar[i] + '[' + std::to_string(refidx[i]) + "] ";
  auto name = kw::amr_refvar::name() + " & id(s)";
  name[0] = static_cast< char >( std::toupper( name[0] ) );
  item( name, c );
}

void InciterPrint::edgeref( const std::vector< std::size_t >& edgenodes )
// *****************************************************************************
// Print initial mesh refinement edge-node pairs
// *****************************************************************************
{
   if (edgenodes.empty()) return;

   std::string c;
   for (auto i : edgenodes) c += std::to_string(i) + ' ';
   auto name = kw::amr_edgelist::name();
   name[0] = static_cast< char >( std::toupper( name[0] ) );
   item( name, c );
}

void InciterPrint::eqlegend()
// *****************************************************************************
// Print PDE factory legend
// *****************************************************************************
{
  section( "PDE factory legend, policy codes" );

  static_assert( tk::HasTypedef_code_v< kw::physics::info >,
                 "Policy code undefined for keyword" );
  static_assert( tk::HasTypedef_code_v< kw::problem::info >,
                 "Policy code undefined for keyword" );

  raw( m_item_indent + kw::physics::name() +
       " (policy code: " + *kw::physics::code() + ")\n" );
  brigand::for_each< ctr::Physics::keywords >( echoPolicies( this ) );

  raw( m_item_indent + kw::problem::name() +
       " (policy code: " + *kw::problem::code() + ")\n" );
  brigand::for_each< ctr::Problem::keywords >( echoPolicies( this ) );
}
