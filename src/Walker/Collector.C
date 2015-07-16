//******************************************************************************
/*!
  \file      src/Walker/Collector.C
  \author    J. Bakosi
  \date      Sat 11 Jul 2015 12:12:12 PM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Charm++ module interface file for collecting contributions from
             Integrators
  \details   Charm++ module interface file for collecting contributions from
             Integrators.
*/
//******************************************************************************

#include "Collector.h"

using walker::Collector;

void
Collector::chareOrd( const std::vector< tk::real >& ord )
//******************************************************************************
// Chares contribute ordinary moments
//! \param[in] ord Vector of partial sums for the estimation of ordinary moments
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
//******************************************************************************
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
//******************************************************************************
// Chares contribute central moments
//! \param[in] cen Vector of partial sums for the estimation of central moments
//! \note This function does not have to be declared as a Charm++ entry
//!   method since it is always called by chares on the same PE.
//! \author J. Bakosi
//******************************************************************************
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "collector.def.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif
