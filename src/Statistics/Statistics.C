//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Sat 05 Jul 2014 08:53:39 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************

#include <cstring>
#include <algorithm>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Base.h>
#include <Statistics.h>

using quinoa::Statistics;

Statistics::Statistics( const Base& base, const ParProps& particles ) :
  m_particles( particles ),
  m_nthreads( base.paradigm.ompNthreads() ),
  m_npar( base.control.get< tag::incpar, tag::npar >() ),
  m_nord(0),
  m_ncen(0)
//******************************************************************************
//  Constructor
//! \param[in]  base       Essentials
//! \param[in]  particles  Particle properties
//! \author  J. Bakosi
//******************************************************************************
{
  // Setup std::map< depvar, offset >
  OffsetMap offset;
  boost::mpl::for_each< tags >( depvar( base, offset ) );

  // Prepare for computing moments
  setupOrdinary( base, offset );
  setupCentral( base, offset );
}

void
Statistics::setupOrdinary( const Base& base, const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  // Prepare for computing ordinary moments
  for (const auto& product : base.control.get< tag::stat >()) {
    if (ordinary(product)) {

      m_instOrd.push_back( std::vector< const tk::real* >() );
      m_plotOrdinary.push_back( false );
      m_nameOrdinary.push_back( std::string() );
      m_ordFieldVar.push_back( ctr::FieldVar() );

      for (const auto& term : product) {
        auto it = offset.find( term.var );
        Assert( it != end( offset ), "No such depvar" );
        // Put in starting address of instantaneous variable
        m_instOrd[ m_nord ].push_back(
          m_particles.cptr( term.field, it->second ) );
        if (term.plot) m_plotOrdinary.back() = true;
        // Put in term name+field
        m_nameOrdinary.back() +=
          m_ordFieldVar.back() = ctr::FieldVar( term.var, term.field );
      }

      ++m_nord;
    }
  }
}

void
Statistics::setupCentral( const Base& base, const OffsetMap& offset )
//******************************************************************************
//  Prepare for computing central moments
//! \author J. Bakosi
//******************************************************************************
{
  if (m_nord) {
    // Storage for all the required ordinary moments
    // +1 for each thread's 0 as center for ordinary moments
    m_ordinary = tk::make_unique< tk::real[] >( m_nthreads*(m_nord+1) );

    // Put in zero as index of center for ordinary moments in central products
    m_ordinary[m_nord] = 0.0;

    // Prepare for computing central moments
    for (const auto& product : base.control.get< tag::stat >()) {
      if (!ordinary(product)) {

        m_instCen.push_back( std::vector< const tk::real* >() );
        m_center.push_back( std::vector< const tk::real* >() );
        m_nameCentral.push_back( std::string() );

        for (const auto& term : product) {
          auto it = offset.find( term.var );
          Assert( it != end( offset ), "No such depvar" );
          // Put in starting address of instantaneous variable
          m_instCen[m_ncen].push_back(
            m_particles.cptr( term.field, it->second ) );
          // Put in index of center for central, m_nord for ordinary moment
          m_center[m_ncen].push_back(
            m_ordinary.get() + (!isupper(term.var) ? mean(term) : m_nord));
          m_nameCentral.back() += ctr::FieldVar( term.var, term.field );
        }

        ++m_ncen;
      }
    }

    if (m_ncen) {
      // Storage for all the required central moments
      m_central = tk::make_unique< tk::real[] >( m_nthreads*m_ncen );
    }
  }
}

bool
Statistics::ordinary( const std::vector< ctr::Term >& product ) const
//******************************************************************************
//  Find out whether product only contains ordinary moment terms
//! \param[in]  product   Vector of terms
//! \author J. Bakosi
//******************************************************************************
{
  // If and only if all terms are ordinary, the product is ordinary
  bool ord = true;
  for (auto& term : product) {
    if (term.moment == ctr::Moment::CENTRAL)
      ord = false;
  }
  return ord;
}

int
Statistics::mean( const ctr::Term& term ) const
//******************************************************************************
//  Return mean for fluctuation
//! \param[in]  term      Term (a fluctuation) whose mean to search for
//! \author J. Bakosi
//******************************************************************************
{
  auto size = m_ordFieldVar.size();
  for (decltype(size) i=0; i<size; ++i) {
    if (m_ordFieldVar[i].var == toupper(term.var) &&
        m_ordFieldVar[i].field == term.field) {
       return i;
    }
  }

  Throw( std::string("Cannot find mean for variable ") + term );
}

bool
Statistics::plotOrdinary( int m ) const
//******************************************************************************
//  Find out whether ordinary moment is to be plotted
//! \param[in]  m         Moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m < m_nord, "Request for unavailable ordinary moment" );
  return m_plotOrdinary[ m ];
}

const std::string&
Statistics::nameOrdinary( int m ) const
//******************************************************************************
//  Return the name of ordinary moment
//! \param[in]  m         Ordinary-moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m < m_nord, "Request for unavailable ordinary moment" );
  return m_nameOrdinary[ m ];
}

const std::string&
Statistics::nameCentral( int m ) const
//******************************************************************************
//  Return the name of central moment
//! \param[in]  m         Central-moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert( m < m_ncen, "Request for unavailable central moment" );
  return m_nameCentral[ m ];
}

void
Statistics::estimateOrdinary()
//******************************************************************************
//  Estimate ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    auto tid = omp_get_thread_num();
    #else
    auto tid = 0;
    #endif

    // Zero ordinary moment accumulators
    memset( m_ordinary.get() + tid*(m_nord+1), 0, m_nord*sizeof(tk::real) );

    // Accumulate ordinary moments
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (uint64_t p=0; p<m_npar; ++p) {
      for (int i=0; i<m_nord; ++i) {
        auto prod = m_particles.cvar( m_instOrd[i][0], p );
        auto s = m_instOrd[i].size();
        for (decltype(s) j=1; j<s; ++j) {
          prod *= m_particles.cvar( m_instOrd[i][j], p );
        }
        m_ordinary[tid*(m_nord+1) + i] += prod;
      }
    }
  }

  // Collect ordinary moments from all threads
  for (uint64_t p=1; p<m_nthreads; ++p) {
    for (int i=0; i<m_nord; ++i) {
      m_ordinary[i] += m_ordinary[p*(m_nord+1) + i];
    }
  }

  // Finish computing ordinary moments
  for (int i=0; i<m_nord; ++i) {
    m_ordinary[i] /= m_npar;
  }
}

void
Statistics::estimateCentral()
//******************************************************************************
//  Estimate central moments
//! \author J. Bakosi
//******************************************************************************
{
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    auto tid = omp_get_thread_num();
    #else
    auto tid = 0;
    #endif

    // Zero central moment accumulators
    memset(m_central.get() + tid*m_ncen, 0, m_ncen*sizeof(tk::real));

    // Accumulate central moments
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (uint64_t p=0; p<m_npar; ++p) {
      for (int i=0; i<m_ncen; ++i) {
        auto prod = m_particles.cvar( m_instCen[i][0], p );
        auto s = m_instCen[i].size();
        for (decltype(s) j=1; j<s; ++j) {
          prod *= m_particles.cvar( m_instCen[i][j], p ) - *(m_center[i][j]);
        }
        m_central[tid*m_ncen + i] += prod;
      }
    }
  }

  // Collect central moments from all threads
  for (uint64_t p=1; p<m_nthreads; ++p) {
    for (int i=0; i<m_ncen; ++i) {
      m_central[i] += m_central[p*m_ncen + i];
    }
  }

  // Finish computing central moments
  for (int i=0; i<m_ncen; ++i) {
    m_central[i] /= m_npar;
  }
}

void
Statistics::accumulate()
//******************************************************************************
//  Acumulate statistics
//! \author J. Bakosi
//******************************************************************************
{
  if (m_nord) {
    estimateOrdinary();                 // Estimate ordinary moments
    if (m_ncen) estimateCentral();      // Estimate central moments
  }
}
