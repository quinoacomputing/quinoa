// *****************************************************************************
/*!
  \file      src/Walker/Collector.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Charm++ module interface file for collecting contributions from
             Integrators
  \details   Charm++ module interface file for collecting contributions from
             Integrators.
*/
// *****************************************************************************

#include "Collector.h"

namespace walker {

//! \brief Charm++ PDF merger reducer
//! \details This variable is defined here in the .C file and declared as extern
//!   in Collector.h. If instead one defines it in the header (as static),
//!   a new version of the variable is created any time the header file is
//!   included, yielding no compilation nor linking errors. However, that leads
//!   to runtime errors, since Collector::registerPDFMerger(), a Charm++
//!   "initnode" entry method, *may* fill one while contribute() may use the
//!   other (unregistered) one. Result: undefined behavior, segfault, and
//!   formatting the internet ...
CkReduction::reducerType PDFMerger;

}

using walker::Collector;

void
Collector::chareOrd( const std::vector< tk::real >& ord )
// *****************************************************************************
// Chares contribute ordinary moments
//! \param[in] ord Vector of partial sums for the estimation of ordinary moments
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  ++m_nord;

  for (std::size_t i=0; i<m_ordinary.size(); ++i) m_ordinary[i] += ord[i];

  // If all chares on my PE have contributed, send partial sums to host
  if (m_nord == m_nchare) {

    // Create Charm++ callback function for reduction
    CkCallback cb( CkReductionTarget( Distributor, estimateOrd ), m_hostproxy );

    // Contribute partial sums to host via Charm++ reduction
    contribute( static_cast< int >( m_ordinary.size() * sizeof(tk::real) ),
                m_ordinary.data(), CkReduction::sum_double, cb );

    // Zero counters for next collection operation
    m_nord = 0;
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );
  }
}

void
Collector::chareCen( const std::vector< tk::real >& cen )
// *****************************************************************************
// Chares contribute central moments
//! \param[in] cen Vector of partial sums for the estimation of central moments
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  ++m_ncen;

  for (std::size_t i=0; i<m_central.size(); ++i) m_central[i] += cen[i];

  // If all chares on my PE have contributed, send partial sums to host
  if (m_ncen == m_nchare) {

    // Create Charm++ callback function for reduction
    CkCallback cb( CkReductionTarget( Distributor, estimateCen ), m_hostproxy );

    // Contribute partial sums to host via Charm++ reduction
    contribute( static_cast< int >( m_central.size() * sizeof(tk::real) ),
                m_central.data(), CkReduction::sum_double, cb );

    // Zero counters for next collection operation
    m_ncen = 0;
    std::fill( begin(m_central), end(m_central), 0.0 );
  }
}

void
Collector::chareOrdPDF( const std::vector< tk::UniPDF >& updf,
                        const std::vector< tk::BiPDF >& bpdf,
                        const std::vector< tk::TriPDF >& tpdf )
// *****************************************************************************
// Chares contribute ordinary PDFs
//! \param[in] updf Vector of partial sums for the estimation of univariate
//!   ordinary PDFs
//! \param[in] bpdf Vector of partial sums for the estimation of bivariate
//!   ordinary PDFs
//! \param[in] tpdf Vector of partial sums for the estimation of trivariate
//!   ordinary PDFs
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  ++m_nopdf;

  // Add contribution from worker chares to partial sums on my PE
  std::size_t i = 0;
  for (const auto& p : updf) m_ordupdf[i++].addPDF( p );
  i = 0;
  for (const auto& p : bpdf) m_ordbpdf[i++].addPDF( p );
  i = 0;
  for (const auto& p : tpdf) m_ordtpdf[i++].addPDF( p );

  // If all chares on my PE have contributed, send partial sums to host
  if (m_nopdf == m_nchare) {

    // Serialize vector of PDFs to raw stream
    auto stream = tk::serialize( m_ordupdf, m_ordbpdf, m_ordtpdf );

    // Create Charm++ callback function for reduction.
    // Distributor::estimateOrdPDF() will be the final target of the reduction
    // where the results of the reduction will appear.
    CkCallback cb( CkIndex_Distributor::estimateOrdPDF(nullptr), m_hostproxy );

    // Contribute serialized PDFs of partial sums to host via Charm++ reduction
    contribute( stream.first, stream.second.get(), PDFMerger, cb );

    // Zero counters for next collection operation
    m_nopdf = 0;
    for (auto& p : m_ordupdf) p.zero();
    for (auto& p : m_ordbpdf) p.zero();
    for (auto& p : m_ordtpdf) p.zero();
  }
}

void
Collector::chareCenPDF( const std::vector< tk::UniPDF >& updf,
                        const std::vector< tk::BiPDF >& bpdf,
                        const std::vector< tk::TriPDF >& tpdf )
// *****************************************************************************
// Chares contribute central PDFs
//! \param[in] updf Vector of partial sums for the estimation of univariate
//!   central PDFs
//! \param[in] bpdf Vector of partial sums for the estimation of bivariate
//!   central PDFs
//! \param[in] tpdf Vector of partial sums for the estimation of trivariate
//!   central PDFs
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
// *****************************************************************************
{
  ++m_ncpdf;

  // Add contribution from worker chares to partial sums on my PE
  std::size_t i = 0;
  for (const auto& p : updf) m_cenupdf[i++].addPDF( p );
  i = 0;
  for (const auto& p : bpdf) m_cenbpdf[i++].addPDF( p );
  i = 0;
  for (const auto& p : tpdf) m_centpdf[i++].addPDF( p );

  // If all chares on my PE have contributed, send partial sums to host
  if (m_ncpdf == m_nchare) {

    // Serialize vector of PDFs to raw stream
    auto stream = tk::serialize( m_cenupdf, m_cenbpdf, m_centpdf );

    // Create Charm++ callback function for reduction.
    // Distributor::estimateCenPDF() will be the final target of the reduction
    // where the results of the reduction will appear.
    CkCallback cb( CkIndex_Distributor::estimateCenPDF(nullptr), m_hostproxy );

    // Contribute serialized PDFs of partial sums to host via Charm++ reduction
    contribute( stream.first, stream.second.get(), PDFMerger, cb );

    // Zero counters for next collection operation
    m_ncpdf = 0;
    for (auto& p : m_cenupdf) p.zero();
    for (auto& p : m_cenbpdf) p.zero();
    for (auto& p : m_centpdf) p.zero();
  }
}

#include "NoWarning/collector.def.h"
