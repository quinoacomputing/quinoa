//******************************************************************************
/*!
  \file      src/Main/WalkerPrint.C
  \author    J. Bakosi
  \date      Thu 19 Mar 2015 12:04:00 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Walker-specific pretty printer functionality
  \details   Walker-specific pretty printer functionality.
*/
//******************************************************************************

#include <WalkerPrint.h>
#include <Walker/Options/DiffEq.h>

using walker::WalkerPrint;

namespace walker {

extern std::vector< DiffEq > g_diffeqs;

} // walker::

void
WalkerPrint::inthead( const std::string& title,
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
WalkerPrint::statistics( const std::string& title ) const
//******************************************************************************
//  Print statistics and PDFs
//! \param[in] title Section title
//! \author J. Bakosi
//******************************************************************************
{
  if ( !g_inputdeck.get< tag::stat >().empty() ||
       !g_inputdeck.get< tag::pdf >().empty() )
  {
    section( title );
    stats( "Estimated statistical moments" );
    pdfs( "PDFs", tk::ctr::pdf );
  }
}

void
WalkerPrint::diffeqs( const std::string& title, const std::vector< std::vector<
  std::pair< std::string, std::string > > >& info ) const
//******************************************************************************
//  Print configuration of a stack of differential equations
//! \param[in] title Title to use
//! \param[in] info Info vector to use
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

void
WalkerPrint::stats( const std::string& msg ) const
//******************************************************************************
//  Echo statistics container contents if differs from default
//! \param[in] msg Message to print
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !msg.empty(), "Empty message size in WalkerPrint::stats()." );
  const auto& c = g_inputdeck.get< tag::stat >();

  if (!c.empty() && c != g_inputdeck_defaults.get< tag::stat >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : c) m_stream << v;
    m_stream << '\n';
  }
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
//******************************************************************************
//  Echo pdfs container contents if differs from default applying op.
//! \param[in] msg Message to print
//! \param[in] op Functor to use
//! \details See src/Control/StatCtr.h for the definition of functions that may
//!    be passed in as op. Currently, the only example is tk::ctr::pdf.
//! \author J. Bakosi
//******************************************************************************
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

  // Oputput options and settings affecting PDF output
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
    item( fl.group(),
          fl.name( g_inputdeck.get< tag::selected, tag::float_format >() ) );
    item( "Text precision in digits",
          g_inputdeck.get< tag::discr, tag::precision >() );
  }
}
