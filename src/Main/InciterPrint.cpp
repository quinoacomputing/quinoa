// *****************************************************************************
/*!
  \file      src/Main/InciterPrint.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter-specific pretty printer functionality
  \details   Inciter-specific pretty printer functionality.
*/
// *****************************************************************************

#include <regex>

#include <brigand/algorithms/for_each.hpp>

#include "CGPDE.hpp"
#include "InciterPrint.hpp"
#include "Transfer.hpp"

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

void
InciterPrint::couple( const std::vector< Transfer >& transfer,
                      const std::vector< char >& depvar ) const
// *****************************************************************************
//  Print out info on solver coupling
//! \param[in] transfer List of solution transfer steps, describing coupling
//! \param[in] depvar List of dependent variables assigned to solvers
// *****************************************************************************
{
  if (!transfer.empty()) {

    endsubsection();
    subsection( "Solver coupling" );
    std::stringstream steps;
    std::map< char, std::vector< char > > src, dst;

    for (const auto& t : transfer) {
      auto sd = depvar[t.src];
      auto dd = depvar[t.dst];
      steps << sd << '>' << dd << ' ';
      src[ sd ].push_back( dd );
      dst[ dd ].push_back( sd );
    }

    item( "Transfer steps (" + std::to_string(transfer.size()) + ')',
          steps.str() );

    for (const auto& [s,m] : src) {
      std::stringstream name;
      name << "Solver " << s << " is source to";
      item( name.str(), tk::parameters(m) );
    }

    for (const auto& [d,m] : dst) {
      std::stringstream name;
      name << "Solver " << d << " is destination to";
      item( name.str(), tk::parameters(m) );
    }

  }
}

void InciterPrint::refvar( const std::vector< std::string >& rvar,
                           const std::vector< std::size_t >& refidx ) const
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

void InciterPrint::edgeref( const std::vector< std::size_t >& edgenodes ) const
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
