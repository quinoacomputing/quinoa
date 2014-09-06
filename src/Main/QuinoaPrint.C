//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.C
  \author    J. Bakosi
  \date      Thu 04 Sep 2014 08:51:14 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaPrint
  \details   QuinoaPrint
*/
//******************************************************************************

#include <QuinoaPrint.h>
#include <Quinoa/Options/DiffEq.h>

using quinoa::QuinoaPrint;

namespace quinoa {

extern std::vector< DiffEq > g_diffeqs;

} // quinoa::

void
QuinoaPrint::stats( const std::string& msg, std::function< std::ostream& (
  std::ostream&, const std::vector< ctr::Term >& ) > op ) const
//******************************************************************************
//  Echo requested statistics if differs from default
//! \author J. Bakosi
//******************************************************************************
{
  //! Fields of std::vector< std::vector< struct{field, name, plot} > > must
  //! exist. See src/Control/Quinoa/InputDeck/Types.h for the definition of
  //! operators, << and <<=, aliased with 'estimated' and 'requested',
  //! respectively, for outputing estimated and requested std::vector< Term >.
  if (g_inputdeck.get<tag::stat>() != g_inputdeck_defaults.get<tag::stat>()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for ( auto& v : g_inputdeck.get< tag::stat >() ) op( m_stream, v );
    m_stream << '\n';
  }
}

void
QuinoaPrint::diffeqs( const std::string& title,
  const std::vector< std::vector< std::pair< std::string, std::string > > >&
    info ) const
//******************************************************************************
//  Print configuration of a stack of differential equations
//! \author J. Bakosi
//******************************************************************************
{
  if ( !info.empty() ) {
    std::stringstream ss;
    ss << title << " (" << g_diffeqs.size() << ")";
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
