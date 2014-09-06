//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Fri 05 Sep 2014 12:35:30 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Computing ordinary and central moments
*/
//******************************************************************************

#include <Statistics.h>

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

  // Use the above offset to prepare for computing all requested ordinary and
  // central moments
  setupOrdinary( offset );
  setupCentral( offset );
}

void
Statistics::setupOrdinary( const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  for (const auto& product : g_inputdeck.get< tag::stat >()) {
    if (ordinary( product )) {

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
}

void
Statistics::setupCentral( const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing central moments
//! \author J. Bakosi
//******************************************************************************
{
  // Central moments can only be estimated if there are ordinary moments
  if (m_nord) {
    // Storage for all the required ordinary moments
    // +1 for 0 as center for ordinary moments
    m_ordinary.resize( m_nord + 1 );

    // Put in zero as center for ordinary moments in central products
    m_ordinary[ m_nord ] = 0.0;

    for (const auto& product : g_inputdeck.get< tag::stat >()) {
      if (central( product )) {

        m_instCen.emplace_back( std::vector< const tk::real* >() );
        m_center.emplace_back( std::vector< const tk::real* >() );

        for (const auto& term : product) {
          auto o = offset.find( term.var );
          Assert( o != end( offset ), "No such depvar" );
          // Put in starting address of instantaneous variable
          m_instCen.back().push_back( m_particles.cptr(term.field, o->second) );
          // Put in index of center for central, m_nord for ordinary moment
          m_center.back().push_back(
            &m_ordinary[0] + (std::islower(term.var) ? mean(term) : m_nord) );
        }

        ++m_ncen;
      }
    }

    // Allocate storage for all required central moments
    m_central.resize( m_ncen );
  }
}

int
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
    m_ordinary = ord;

    // Zero central moment accumulators
    std::fill( begin(m_central), end(m_central), 0.0 );

    // Accumulate sum for central moments. This is a partial sum, so no division
    // by the number of samples.
    const auto npar = m_particles.npar();
    for (auto p=decltype(npar){0}; p<npar; ++p) {
      for (int i=0; i<m_ncen; ++i) {
        auto prod = m_particles.cvar( m_instCen[i][0], p ) - *(m_center[i][0]);
        const auto s = m_instCen[i].size();
        for (auto j=decltype(s){1}; j<s; ++j) {
          prod *= m_particles.cvar( m_instCen[i][j], p ) - *(m_center[i][j]);
        }
        m_central[i] += prod;
      }
    }
  }
}
