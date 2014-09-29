//******************************************************************************
/*!
  \file      src/Main/QuinoaPrint.C
  \author    J. Bakosi
  \date      Mon 22 Sep 2014 03:59:58 PM MDT
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
QuinoaPrint::inthead( const std::string& title,
                      const std::string& name,
                      const std::string& legend,
                      const std::string& head ) const
//******************************************************************************
//  Print time integration header
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
//! \author J. Bakosi
//******************************************************************************
{
  if ( !g_inputdeck.get< tag::stat >().empty() ||
       !g_inputdeck.get< tag::pdf >().empty() )
  {
    section( title );
    stats( "Requested statistical moments", ctr::requested );
    stats( "Triggered statistical moments", ctr::triggered );
    stats( "Estimated statistical moments", ctr::estimated );
    pdfs( "PDFs", ctr::pdf );
  }
}

void
QuinoaPrint::diffeqs( const std::string& title, const std::vector< std::vector<
  std::pair< std::string, std::string > > >& info ) const
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

void
QuinoaPrint::stats( const std::string& msg, std::function< std::ostream& (
  std::ostream&, const std::vector< ctr::Term >& ) > op ) const
//******************************************************************************
//  Echo statistics container contents if differs from default applying op.
//! \details See src/Control/Quinoa/InputDeck/Types.h for the definition of
//! functions that may be passed in as op. Examples are 'estimated',
//! 'requested', and 'triggered'. The operation given by the template
//! argument and is a function pointer specifying an stream-output operator
//! for a std::vector< ctr::Term >.
//! \author J. Bakosi
//******************************************************************************
{
  Assert( !msg.empty(), "Empty message size in QuinoaPrint::stats()." );
  const auto& c = g_inputdeck.get< tag::stat >();

  if (!c.empty() && c != g_inputdeck_defaults.get< tag::stat >()) {
    m_stream << m_item_name_fmt % m_item_indent % msg;
    for (auto& v : c) op( m_stream, v );
    m_stream << '\n';
  }
}

void
QuinoaPrint::pdfs( const std::string& msg,
                   std::function<
                     std::ostream& ( std::ostream&,
                                     const std::vector< ctr::Term >&,
                                     const std::vector< tk::real >&,
                                     const std::string&,
                                     const std::vector< tk::real >& ext ) > op )
const
//******************************************************************************
//  Echo pdfs container contents if differs from default applying op.
//! \details See src/Control/Quinoa/InputDeck/Types.h for the definition of
//! functions that may be passed in as op. Examples are 'estimated',
//! 'requested', and 'triggered'. The operation given by the template
//! argument and is a function pointer specifying an stream-output operator
//! for a std::vector< ctr::Term >.
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

  // Print PDF output options
  ctr::PDFFile f;
  item( f.group(),
        f.name( g_inputdeck.get< tag::selected, tag::pdffiletype >() ) );
  ctr::PDFPolicy p;
  item( p.group(),
        p.name( g_inputdeck.get< tag::selected, tag::pdfpolicy >() ) );
  ctr::PDFCentering e;
  item( e.group(),
        e.name( g_inputdeck.get< tag::selected, tag::pdfctr >() ) );
}
