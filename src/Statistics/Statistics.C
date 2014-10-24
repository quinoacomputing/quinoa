//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Sun 19 Oct 2014 09:47:06 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Computing ordinary and central moments
*/
//******************************************************************************

#include <Statistics.h>
#include <flip_map.h>

using quinoa::Statistics;

Statistics::Statistics( const ParProps& particles ) :
  m_particles( particles ),
  m_nord( 0 ),
  m_ncen( 0 )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Setup a map that associates each dependent variable of all requested
  // differential equations to their offset in the particle array
  OffsetMap offset;
  boost::mpl::for_each< ncomps::tags >( depvar( offset ) );

  // Prepare for computing ordinary and central moments, PDFs
  setupOrdinary( offset );
  setupCentral( offset );
  setupPDF( offset );
}

void
Statistics::setupOrdinary( const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  for (const auto& product : g_inputdeck.get< tag::stat >()) {
    if (ordinary(product)) {

      m_instOrd.emplace_back( std::vector< const tk::real* >() );
      m_ordFieldVar.emplace_back( ctr::FieldVar() );

      for (const auto& term : product) {
        auto o = offset.find( term.var );
        Assert( o != end( offset ), "No such depvar" );
        // Put in starting address of instantaneous variable
        m_instOrd.back().push_back( m_particles.cptr( term.field, o->second ) );
        // Put in term name+field
        m_ordFieldVar.back() = ctr::FieldVar( term.var, term.field );
      }

      ++m_nord;
    }
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
Statistics::setupCentral( const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing central moments
//! \author J. Bakosi
//******************************************************************************
{
  // Central moments can only be estimated about ordinary moments
  if (m_nord) {
    for (const auto& product : g_inputdeck.get< tag::stat >()) {
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
            &m_ordinary[0] + (std::islower(term.var) ? mean(term) : m_nord) );
        }

        ++m_ncen;
      }
    }

    // Allocate storage for all required central moments
    m_central.resize( m_ncen );
  }
}

void
Statistics::setupPDF( const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing PDFs
//! \author J. Bakosi
//******************************************************************************
{
  std::size_t i = 0;
  for (const auto& probability : g_inputdeck.get< tag::pdf >()) {
    if (ordinary(probability)) {

      // Detect number of sample space dimensions and create ordinary PDFs
      const auto& bs = g_inputdeck.get< tag::discr, tag::binsize >()[i++];
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
      const auto& bs = g_inputdeck.get< tag::discr, tag::binsize >()[i++];
      if (bs.size() == 1) {
        m_cenupdf.emplace_back( bs[0] );
        m_instCenUniPDF.emplace_back( std::vector< const tk::real* >() );
        m_ctrUniPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 2) {
        m_cenbpdf.emplace_back( bs );
        m_instCenBiPDF.emplace_back( std::vector< const tk::real* >() );
        m_ctrBiPDF.emplace_back( std::vector< const tk::real* >() );
      } else if (bs.size() == 3) {
        m_ordtpdf.emplace_back( bs );
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
          &m_ordinary[0] + (std::islower(term.var) ? mean(term) : m_nord);
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
Statistics::mean( const ctr::Term& term ) const
//******************************************************************************
//  Return mean for fluctuation
//! \param[in]  term      Term (a fluctuation) whose mean to search for
//! \author J. Bakosi
//******************************************************************************
{
  const auto size = m_ordFieldVar.size();
  for (auto i=decltype(size){0}; i<size; ++i) {
    if (m_ordFieldVar[i].var == std::toupper(term.var) &&
        m_ordFieldVar[i].field == term.field) {
      return i;
    }
  }

  Throw( std::string("Cannot find mean for variable ") + term );
}

void
Statistics::accumulateOrd()
//******************************************************************************
//  Accumulate (i.e., only do the sum for) ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  if (m_nord) {
    // Zero ordinary moment accumulators
    std::fill( begin(m_ordinary), end(m_ordinary), 0.0 );

    // Accumulate sum for ordinary moments. This is a partial sum, so no division
    // by the number of samples.
    const auto npar = m_particles.npar();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      for (int i=0; i<m_nord; ++i) {
       auto prod = m_particles.cvar( m_instOrd[i][0], p );
        const auto s = m_instOrd[i].size();
        for (auto j=decltype(s){1}; j<s; ++j) {
          prod *= m_particles.cvar( m_instOrd[i][j], p );
        }
        m_ordinary[i] += prod;
      }
    }
  }
}

void
Statistics::accumulateCen( const std::vector< tk::real >& ord )
//******************************************************************************
//  Accumulate (i.e., only do the sum for) central moments
//  \param[in]  ord  Ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  if (m_ncen) {
    // Overwrite ordinary moments by those computed across all PEs
    for (std::size_t i=0; i<ord.size(); ++i) m_ordinary[i] = ord[i];

    // Zero central moment accumulators
    std::fill( begin(m_central), end(m_central), 0.0 );

    // Accumulate sum for central moments. This is a partial sum, so no division
    // by the number of samples.
    const auto npar = m_particles.npar();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      for (int i=0; i<m_ncen; ++i) {
        auto prod = m_particles.cvar( m_instCen[i][0], p ) - *(m_ctr[i][0]);
        const auto s = m_instCen[i].size();
        for (auto j=decltype(s){1}; j<s; ++j) {
          prod *= m_particles.cvar( m_instCen[i][j], p ) - *(m_ctr[i][j]);
        }
        m_central[i] += prod;
      }
    }
  }
}

void
Statistics::accumulateOrdPDF()
//******************************************************************************
//  Accumulate (i.e., only do the sum for) ordinary PDFs
//! \author J. Bakosi
//******************************************************************************
{
  if (!m_ordupdf.empty() || !m_ordbpdf.empty() || !m_ordtpdf.empty()) {
    // Zero PDF accumulators
    for (auto& pdf : m_ordupdf) pdf.zero();
    for (auto& pdf : m_ordbpdf) pdf.zero();
    for (auto& pdf : m_ordtpdf) pdf.zero();

    // Accumulate partial sum for PDFs
    const auto npar = m_particles.npar();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      std::size_t i = 0;
      // Accumulate partial sum for univariate PDFs
      for (auto& pdf : m_ordupdf) {
        pdf.add( m_particles.cvar( m_instOrdUniPDF[i++][0], p ) );
      }
      // Accumulate partial sum for bivariate PDFs
      i = 0;
      for (auto& pdf : m_ordbpdf) {
        const auto inst = m_instOrdBiPDF[i++];
        pdf.add( {{ m_particles.cvar( inst[0], p ),
                    m_particles.cvar( inst[1], p ) }} );
      }
      // Accumulate partial sum for trivariate PDFs
      i = 0;
      for (auto& pdf : m_ordtpdf) {
        const auto inst = m_instOrdTriPDF[i++];
        pdf.add( {{ m_particles.cvar( inst[0], p ),
                    m_particles.cvar( inst[1], p ),
                    m_particles.cvar( inst[2], p ) }} );
      }
    }
  }
}

void
Statistics::accumulateCenPDF( const std::vector< tk::real >& ord )
//******************************************************************************
//  Accumulate (i.e., only do the sum for) central PDFs
//! \author J. Bakosi
//******************************************************************************
{
  if (!m_cenupdf.empty() || !m_cenbpdf.empty() || !m_centpdf.empty()) {
    // Overwrite ordinary moments by those computed across all PEs
    for (std::size_t i=0; i<ord.size(); ++i) m_ordinary[i] = ord[i];

    // Zero PDF accumulators
    for (auto& pdf : m_cenupdf) pdf.zero();
    for (auto& pdf : m_cenbpdf) pdf.zero();
    for (auto& pdf : m_centpdf) pdf.zero();

    // Accumulate partial sum for PDFs
    const auto npar = m_particles.npar();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      std::size_t i = 0;
      // Accumulate partial sum for univariate PDFs
      for (auto& pdf : m_cenupdf) {
        pdf.add(
          m_particles.cvar( m_instCenUniPDF[i][0], p ) - *(m_ctrUniPDF[i][0]) );
        ++i;
      }
      // Accumulate partial sum for bivariate PDFs
      i = 0;
      for (auto& pdf : m_cenbpdf) {
        const auto& inst = m_instCenBiPDF[i];
        const auto& cen = m_ctrBiPDF[i];
        pdf.add( {{ m_particles.cvar( inst[0], p ) - *(cen[0]),
                    m_particles.cvar( inst[1], p ) - *(cen[1]) }} );
        ++i;
      }
      // Accumulate partial sum for trivariate PDFs
      i = 0;
      for (auto& pdf : m_centpdf) {
        const auto inst = m_instCenTriPDF[i];
        const auto& cen = m_ctrTriPDF[i];
        pdf.add( {{ m_particles.cvar( inst[0], p ) - *(cen[0]),
                    m_particles.cvar( inst[1], p ) - *(cen[1]),
                    m_particles.cvar( inst[2], p ) - *(cen[2]) }} );
        ++i;
      }
    }
  }
}
