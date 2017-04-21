// *****************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Statistics class definition
  \details   This file implements a statistics class that can be used to
    estimate statistics from an ensemble. Supported at this time are ordinary
    and central statistical moments of arbitrary-length products and arbitrary
    number of 1D, 2D, and 3D probability density functions (PDF) with sample
    spaces of ordinary and/or central sample space variables. See the header
    file documentation for more information on the nomenclature.
*/
// *****************************************************************************

#include <map>
#include <iterator>
#include <utility>
#include <algorithm>
#include <iosfwd>
#include <cctype>

#include "Types.h"
#include "Exception.h"
#include "Statistics.h"
#include "Particles.h"
#include "SystemComponents.h"
#include "UniPDF.h"
#include "BiPDF.h"
#include "TriPDF.h"

using tk::Statistics;

Statistics::Statistics( const tk::Particles& particles,
                        const ctr::OffsetMap& offset,
                        const std::vector< ctr::Product >& stat,
                        const std::vector< ctr::Probability >& pdf,
                        const std::vector< std::vector< tk::real > >& binsize )
  : m_particles( particles ),
    m_instOrd(),
    m_ordinary(),
    m_ordTerm(),
    m_nord( 0 ),
    m_instCen(),
    m_central(),
    m_ctr(),
    m_ncen( 0 ),
    m_instOrdUniPDF(),
    m_ordupdf(),
    m_instCenUniPDF(),
    m_cenupdf(),
    m_ctrUniPDF(),
    m_instOrdBiPDF(),
    m_ordbpdf(),
    m_instCenBiPDF(),
    m_cenbpdf(),
    m_ctrBiPDF(),
    m_instOrdTriPDF(),
    m_ordtpdf(),
    m_instCenTriPDF(),
    m_centpdf(),
    m_ctrTriPDF()
// *****************************************************************************
//  Constructor
//! \param[in] particles Particles data to estimate from
//! \param[in] offset Map of offsets in memory to address variable fields
//! \param[in] stat List of requested statistical moments
//! \param[in] pdf List of requested probability density functions (PDF)
//! \param[in] binsize List of binsize vectors configuring the PDF estimators
//! \author  J. Bakosi
// *****************************************************************************
{
  // Prepare for computing ordinary and central moments, PDFs
  setupOrdinary( offset, stat );
  setupCentral( offset, stat );
  setupPDF( offset, pdf, binsize );
}

void
Statistics::setupOrdinary( const ctr::OffsetMap& offset,
                           const std::vector< ctr::Product >& stat )
// *****************************************************************************
//  Prepare for computing ordinary moments
//! \param[in] offset Map of offsets in memory to address variable fields
//! \param[in] stat List of requested statistical moments
//! \author J. Bakosi
// *****************************************************************************
{
  for (const auto& product : stat)
    if (ordinary(product)) {

      m_instOrd.emplace_back( std::vector< const tk::real* >() );

      int i = 0;
      for (const auto& term : product) {
        auto o = offset.find( term.var );
        Assert( o != end( offset ), "No such depvar" );
        // Put in starting address of instantaneous variable
        m_instOrd.back().push_back( m_particles.cptr( term.field, o->second ) );
        // Collect all means of estimated statistics in a linear vector; this
        // will be used to find means for fluctuations. Thus only collect single
        // terms, i.e., <Y1>, <Y2>, etc., but not <Y1Y2>, etc.
        if (i==0) m_ordTerm.push_back( term );
        ++i;
      }

      // Increase number of ordinary moments by one
      m_ordinary.push_back( 0.0 );
      // Count up orindary moments
      ++m_nord;
    }

  // Put in a zero as the last ordinary moment. This will be used as the center
  // about which central moments are computed. If this is not needed, e.g.,
  // because there is no central moments or not central PDFs are requested, this
  // is small, unused, and harmless.
  if (m_nord) {
    // Storage for all the required ordinary moments
    // +1 for 0 as center for ordinary moments in computing central moments
    m_ordinary.resize( m_nord + 1 );
    // Put in zero as center for ordinary moments in central products
    m_ordinary[ m_nord ] = 0.0;
  }
}

void
Statistics::setupCentral( const ctr::OffsetMap& offset,
                          const std::vector< ctr::Product >& stat )
// *****************************************************************************
//  Prepare for computing central moments
//! \param[in] offset Map of offsets in memory to address variable fields
//! \param[in] stat List of requested statistical moments
//! \author J. Bakosi
// *****************************************************************************
{
  // Central moments can only be estimated about ordinary moments
  if (m_nord)
    for (const auto& product : stat) {
      if (central(product)) {

        m_instCen.emplace_back( std::vector< const tk::real* >() );
        m_ctr.emplace_back( std::vector< const tk::real* >() );

        for (const auto& term : product) {
          auto o = offset.find( term.var );
          Assert( o != end( offset ), "No such depvar" );
          // Put in starting address of instantaneous variable
          m_instCen.back().push_back( m_particles.cptr(term.field, o->second) );
          // Put in index of center for central, m_nord for ordinary moment
          m_ctr.back().push_back(
           m_ordinary.data() + (std::islower(term.var) ? mean(term) : m_nord) );
        }

        // Increase number of central moments by one
        m_central.push_back( 0.0 );
        // Count up central moments
        ++m_ncen;
      }
    }
}

void
Statistics::setupPDF( const ctr::OffsetMap& offset,
                      const std::vector< ctr::Probability >& pdf,
                      const std::vector< std::vector< tk::real > >& binsize )
// *****************************************************************************
//  Prepare for computing PDFs
//! \param[in] offset Map of offsets in memory to address variable fields
//! \param[in] pdf List of requested probability density functions (PDF)
//! \param[in] binsize List of binsize vectors configuring the PDF estimators
//! \author J. Bakosi
// *****************************************************************************
{
  std::size_t i = 0;
  for (const auto& probability : pdf) {
    if (ordinary(probability)) {

      // Detect number of sample space dimensions and create ordinary PDFs
      const auto& bs = binsize[i++];
      if (bs.size() == 1) {
        m_ordupdf.emplace_back( bs[0] );
        m_instOrdUniPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 2) {
        m_ordbpdf.emplace_back( bs );
        m_instOrdBiPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 3) {
        m_ordtpdf.emplace_back( bs );
        m_instOrdTriPDF.emplace_back( std::vector< const tk::real* >() );
      }

      // Put in starting addresses of instantaneous variables
      for (const auto& term : probability) {
        auto o = offset.find( term.var );
        Assert( o != end( offset ), "No such depvar" );
        const tk::real* iptr = m_particles.cptr( term.field, o->second );
        if (bs.size() == 1) m_instOrdUniPDF.back().push_back( iptr );
        else if (bs.size() == 2) m_instOrdBiPDF.back().push_back( iptr );
        else if (bs.size() == 3) m_instOrdTriPDF.back().push_back( iptr );
      }

    } else { // if central PDF

      // Detect number of sample space dimensions and create central PDFs,
      // create new storage for instantaneous variable pointer, create new
      // storage for center pointer
      const auto& bs = binsize[i++];
      if (bs.size() == 1) {
        m_cenupdf.emplace_back( bs[0] );
        m_instCenUniPDF.emplace_back( std::vector< const tk::real* >() );
        m_ctrUniPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 2) {
        m_cenbpdf.emplace_back( bs );
        m_instCenBiPDF.emplace_back( std::vector< const tk::real* >() );
        m_ctrBiPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 3) {
        m_centpdf.emplace_back( bs );
        m_instCenTriPDF.emplace_back( std::vector< const tk::real* >() );
        m_ctrTriPDF.emplace_back( std::vector< const tk::real* >() );
      }

      // Put in starting address of instantaneous variables
      for (const auto& term : probability) {
        auto o = offset.find( term.var );
        Assert( o != end( offset ), "No such depvar" );
        // Put in starting address of instantaneous variable as well as index
        // of center for central, m_nord for ordinary moment
        const tk::real* iptr = m_particles.cptr( term.field, o->second );
        const tk::real* cptr =
          m_ordinary.data() + (std::islower(term.var) ? mean(term) : m_nord);
        if (bs.size() == 1) {
          m_instCenUniPDF.back().push_back( iptr );
          m_ctrUniPDF.back().push_back( cptr );
        } else if (bs.size() == 2) {
          m_instCenBiPDF.back().push_back( iptr );
          m_ctrBiPDF.back().push_back( cptr );
        } else if (bs.size() == 3) {
          m_instCenTriPDF.back().push_back( iptr );
          m_ctrTriPDF.back().push_back( cptr );
        }
      }
    }
  }
}

std::size_t
Statistics::mean( const tk::ctr::Term& term ) const
// *****************************************************************************
//  Return mean for fluctuation
//! \param[in] term Term (a fluctuation) whose mean to search for
//! \return Index to mean
//! \author J. Bakosi
// *****************************************************************************
{
  const auto size = m_ordTerm.size();
  for (auto i=decltype(size){0}; i<size; ++i) {
    if (m_ordTerm[i].var == std::toupper(term.var) &&
        m_ordTerm[i].field == term.field) {
      return i;
    }
  }
  Throw( std::string("Cannot find mean for variable ") + term );
}

void
Statistics::accumulateOrd()
// *****************************************************************************
//  Accumulate (i.e., only do the sum for) ordinary moments
//! \author J. Bakosi
// *****************************************************************************
{
  if (m_nord) {
    // Zero ordinary moment accumulators
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );

    // Accumulate sum for ordinary moments. This is a partial sum, so no
    // division by the number of samples.
    const auto npar = m_particles.nunk();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      for (std::size_t i=0; i<m_nord; ++i) {
       auto prod = m_particles.var( m_instOrd[i][0], p );
        const auto s = m_instOrd[i].size();
        for (auto j=decltype(s){1}; j<s; ++j) {
          prod *= m_particles.var( m_instOrd[i][j], p );
        }
        m_ordinary[i] += prod;
      }
    }
  }
}

void
Statistics::accumulateCen( const std::vector< tk::real >& om )
// *****************************************************************************
//  Accumulate (i.e., only do the sum for) central moments
//! \details The ordinary moments container, m_ordinary, is overwritten here
//!   with the argument om, because each of multiple Statistics class objects
//!   (residing on different PEs) only collect their partial sums when
//!   accumulateOrd() is run. By the time the accumulation of the central
//!   moments is started, the ordinary moments have been collected from all
//!   PEs and thus are the same to be passed here on all PEs. For example
//!   client-code, see walker::Distributor.
//! \param[in] om Ordinary moments
//! \author J. Bakosi
// *****************************************************************************
{
  if (m_ncen) {
    // Overwrite ordinary moments by those computed across all PEs
    for (std::size_t i=0; i<om.size(); ++i) m_ordinary[i] = om[i];

    // Zero central moment accumulators
    std::fill( begin(m_central), end(m_central), 0.0 );

    // Accumulate sum for central moments. This is a partial sum, so no division
    // by the number of samples.
    const auto npar = m_particles.nunk();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      for (std::size_t i=0; i<m_ncen; ++i) {
        auto prod = m_particles.var( m_instCen[i][0], p ) - *(m_ctr[i][0]);
        const auto s = m_instCen[i].size();
        for (auto j=decltype(s){1}; j<s; ++j) {
          prod *= m_particles.var( m_instCen[i][j], p ) - *(m_ctr[i][j]);
        }
        m_central[i] += prod;
      }
    }
  }
}

void
Statistics::accumulateOrdPDF()
// *****************************************************************************
//  Accumulate (i.e., only do the sum for) ordinary PDFs
//! \author J. Bakosi
// *****************************************************************************
{
  if (!m_ordupdf.empty() || !m_ordbpdf.empty() || !m_ordtpdf.empty()) {
    // Zero PDF accumulators
    for (auto& pdf : m_ordupdf) pdf.zero();
    for (auto& pdf : m_ordbpdf) pdf.zero();
    for (auto& pdf : m_ordtpdf) pdf.zero();

    // Accumulate partial sum for PDFs
    const auto npar = m_particles.nunk();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      std::size_t i = 0;
      // Accumulate partial sum for univariate PDFs
      for (auto& pdf : m_ordupdf) {
        pdf.add( m_particles.var( m_instOrdUniPDF[i++][0], p ) );
      }
      // Accumulate partial sum for bivariate PDFs
      i = 0;
      for (auto& pdf : m_ordbpdf) {
        const auto inst = m_instOrdBiPDF[i++];
        pdf.add( {{ m_particles.var( inst[0], p ),
                    m_particles.var( inst[1], p ) }} );
      }
      // Accumulate partial sum for trivariate PDFs
      i = 0;
      for (auto& pdf : m_ordtpdf) {
        const auto inst = m_instOrdTriPDF[i++];
        pdf.add( {{ m_particles.var( inst[0], p ),
                    m_particles.var( inst[1], p ),
                    m_particles.var( inst[2], p ) }} );
      }
    }
  }
}

void
Statistics::accumulateCenPDF( const std::vector< tk::real >& om )
// *****************************************************************************
//  Accumulate (i.e., only do the sum for) central PDFs
//! \details The ordinary moments container, m_ordinary, is overwritten here
//!   with the argument om, because each of multiple Statistics class objects
//!   (residing on different PEs) only collect their partial sums when
//!   accumulateOrd() is run. By the time the accumulation of the central
//!   PDFs is started, the ordinary moments have been collected from all
//!   PEs and thus are the same to be passed here on all PEs. For example
//!   client-code, see walker::Distributor.
//! \param[in] om Ordinary moments
//! \author J. Bakosi
// *****************************************************************************
{
  if (!m_cenupdf.empty() || !m_cenbpdf.empty() || !m_centpdf.empty()) {
    // Overwrite ordinary moments by those computed across all PEs
    for (std::size_t i=0; i<om.size(); ++i) m_ordinary[i] = om[i];

    // Zero PDF accumulators
    for (auto& pdf : m_cenupdf) pdf.zero();
    for (auto& pdf : m_cenbpdf) pdf.zero();
    for (auto& pdf : m_centpdf) pdf.zero();

    // Accumulate partial sum for PDFs
    const auto npar = m_particles.nunk();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      std::size_t i = 0;
      // Accumulate partial sum for univariate PDFs
      for (auto& pdf : m_cenupdf) {
        pdf.add(
          m_particles.var( m_instCenUniPDF[i][0], p ) - *(m_ctrUniPDF[i][0]) );
        ++i;
      }
      // Accumulate partial sum for bivariate PDFs
      i = 0;
      for (auto& pdf : m_cenbpdf) {
        const auto& inst = m_instCenBiPDF[i];
        const auto& cen = m_ctrBiPDF[i];
        pdf.add( {{ m_particles.var( inst[0], p ) - *(cen[0]),
                    m_particles.var( inst[1], p ) - *(cen[1]) }} );
        ++i;
      }
      // Accumulate partial sum for trivariate PDFs
      i = 0;
      for (auto& pdf : m_centpdf) {
        const auto inst = m_instCenTriPDF[i];
        const auto& cen = m_ctrTriPDF[i];
        pdf.add( {{ m_particles.var( inst[0], p ) - *(cen[0]),
                    m_particles.var( inst[1], p ) - *(cen[1]),
                    m_particles.var( inst[2], p ) - *(cen[2]) }} );
        ++i;
      }
    }
  }
}
