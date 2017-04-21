// *****************************************************************************
/*!
  \file      src/Main/WalkerPrint.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Walker-specific pretty printer functionality
  \details   Walker-specific pretty printer functionality.
*/
// *****************************************************************************

#include "NoWarning/replace.h"
#include <boost/iterator/iterator_traits.hpp>

#include "Exception.h"
#include "Tags.h"
#include "StatCtr.h"
#include "WalkerPrint.h"
#include "DiffEq.h"
#include "Options/PDFCentering.h"
#include "Options/PDFFile.h"
#include "Options/PDFPolicy.h"
#include "Options/TxtFloatFormat.h"

using walker::WalkerPrint;

namespace walker {

extern std::vector< DiffEq > g_diffeqs;

} // walker::

void
WalkerPrint::inthead( const std::string& t,
                      const std::string& name,
                      const std::string& legend,
                      const std::string& head ) const
// *****************************************************************************
//  Print time integration header
//! \param[in] t Section title
//! \param[in] name Section name
//! \param[in] legend Legend to print
//! \param[in] head Head to append
//! \author J. Bakosi
// *****************************************************************************
{
  section( t, name );
  std::string l( legend );
  boost::replace_all( l, "\n", "\n" + m_item_indent );
  raw( m_item_indent + l + head );
}

void
WalkerPrint::statistics( const std::string& t ) const
// *****************************************************************************
//  Print statistics and PDFs
//! \param[in] t Section title
//! \author J. Bakosi
// *****************************************************************************
{
  if ( !g_inputdeck.get< tag::stat >().empty() ||
       !g_inputdeck.get< tag::pdf >().empty() )
  {
    section( t );
    stats( "Estimated statistical moments" );
    pdfs( "Estimated PDFs", tk::ctr::pdf );
  }
}

void
WalkerPrint::diffeqs( const std::string& t, const std::vector< std::vector<
  std::pair< std::string, std::string > > >& info ) const
// *****************************************************************************
//  Print configuration of a stack of differential equations
//! \param[in] t Title to use
//! \param[in] info Info vector to use
//! \author J. Bakosi
// *****************************************************************************
{
  if ( !info.empty() ) {
    std::stringstream ss;
    ss << t << " (" << g_diffeqs.size() << ")";
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
WalkerPrint::stats( const std::string& msg ) const
// *****************************************************************************
//  Echo statistics container contents if differs from default
//! \param[in] msg Message to print
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( !msg.empty(), "Empty message size in WalkerPrint::stats()." );
  const auto& c = g_inputdeck.get< tag::stat >();

  if (!c.empty() && c != g_inputdeck_defaults.get< tag::stat >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : c) m_stream << v;
    m_stream << '\n';
  }

  // Output options and settings affecting statistics output
  tk::ctr::TxtFloatFormat fl;
  item( "Stats " + fl.group(),
        fl.name( g_inputdeck.get< tag::flformat, tag::stat >() ) );
  item( "Stats text precision, digits",
        g_inputdeck.get< tag::prec, tag::stat >() );
}

void
WalkerPrint::pdfs( const std::string& msg,
                   std::function<
                     std::ostream& ( std::ostream&,
                                     const std::vector< tk::ctr::Term >&,
                                     const std::vector< tk::real >&,
                                     const std::string&,
                                     const std::vector< tk::real >& ext ) > op )
const
// *****************************************************************************
//  Echo pdfs container contents if differs from default applying op.
//! \param[in] msg Message to print
//! \param[in] op Functor to use
//! \details See src/Control/StatCtr.h for the definition of functions that may
//!    be passed in as op. Currently, the only example is tk::ctr::pdf.
//! \author J. Bakosi
// *****************************************************************************
{
  Assert( !msg.empty(), "Empty message size in WalkerPrint::vec()." );

  const auto& c = g_inputdeck.get< tag::pdf >();
  const auto& b = g_inputdeck.get< tag::discr, tag::binsize >();
  const auto& n = g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >();
  const auto& x = g_inputdeck.get< tag::discr, tag::extent >();

  Assert( (c.size() == b.size()) &&
          (c.size() == n.size()) &&
          (c.size() == x.size()),
          "Number of PDFs, number of binsizes vector, number of PDF names, and "
          "number of extents vector must all equal in WalkerPrint::pdfs()." );

  if (!c.empty() && c != g_inputdeck_defaults.get< tag::pdf >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (std::size_t i=0; i<c.size(); ++i) {
      op( m_stream, c[i], b[i], n[i], x[i] );
    }
    m_stream << '\n';
  }

  // Output options and settings affecting PDF output
  if (!c.empty()) {
    tk::ctr::PDFFile f;
    item( f.group(),
          f.name( g_inputdeck.get< tag::selected, tag::pdffiletype >() ) );
    tk::ctr::PDFPolicy p;
    item( p.group(),
          p.name( g_inputdeck.get< tag::selected, tag::pdfpolicy >() ) );
    tk::ctr::PDFCentering e;
    item( e.group(),
          e.name( g_inputdeck.get< tag::selected, tag::pdfctr >() ) );
    tk::ctr::TxtFloatFormat fl;
    item( "PDF text " + fl.group(),
          fl.name( g_inputdeck.get< tag::flformat, tag::pdf >() ) );
    item( "PDF text precision, digits",
          g_inputdeck.get< tag::prec, tag::pdf >() );
  }
}
