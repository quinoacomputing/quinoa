//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.C
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 12:03:14 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Quinoa-specific pretty printer functionality
  \details   Quinoa-specific pretty printer functionality.
*/
//******************************************************************************

#include <QuinoaPrint.h>

using quinoa::QuinoaPrint;

void
QuinoaPrint::inthead( const std::string& title,
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
QuinoaPrint::statistics( const std::string& title ) const
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
QuinoaPrint::stats( const std::string& msg ) const
//******************************************************************************
//  Echo statistics container contents if differs from default
//! \param[in] msg Message to print
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !msg.empty(), "Empty message size in QuinoaPrint::stats()." );
  const auto& c = g_inputdeck.get< tag::stat >();

  if (!c.empty() && c != g_inputdeck_defaults.get< tag::stat >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : c) m_stream << v;
    m_stream << '\n';
  }
}

void
QuinoaPrint::pdfs( const std::string& msg,
                   std::function<
                     std::ostream& ( std::ostream&,
                                     const std::vector< tk::ctr::Term >&,
                                     const std::vector< tk::real >&,
                                     const std::string&,
                                     const std::vector< tk::real >& ext ) > op )
const
//******************************************************************************
//  Echo pdfs container contents if differs from default applying op().
//! \param[in] msg Message to print
//! \param[in] op Functor to use
//! \details See src/Control/StatCtr.h for the definition of functions that may
//!    be passed in as op. Currently, the only example is tk::ctr::pdf.
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !msg.empty(), "Empty message size in QuinoaPrint::vec()." );

  const auto& c = g_inputdeck.get< tag::pdf >();
  const auto& b = g_inputdeck.get< tag::discr, tag::binsize >();
  const auto& n = g_inputdeck.get< tag::cmd, tag::io, tag::pdfnames >();
  const auto& x = g_inputdeck.get< tag::discr, tag::extent >();

  Assert( (c.size() == b.size()) &&
          (c.size() == n.size()) &&
          (c.size() == x.size()),
          "Number of PDFs, number of binsizes vector, number of PDF names, and "
          "number of extents vector must all equal in QuinoaPrint::pdfs()." );

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
