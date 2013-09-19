//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Thu Sep 19 13:31:18 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
#include <Paradigm.h>
#include <QuinoaControl.h>

using namespace quinoa;

Statistics::Statistics(const Base& base, const real* const particles)
//******************************************************************************
//  Constructor
//! \param[in]  base       Essentials
//! \param[in]  particles  Particles
//! \author  J. Bakosi
//******************************************************************************
try :
  m_base(base),
  m_nthread(base.paradigm.nthread()),
  m_npar(base.control.get<ctr::component, ctr::npar>()),
  m_particles(particles),
  m_nprop(base.control.nprop()),
  m_statistics(base.control.get<ctr::stats>()),
  m_instOrd(),
  m_ordinary(),
  m_ordFieldName(),
  m_nameOrdinary(),
  m_nord(0),
  m_instCen(),
  m_central(),
  m_nameCentral(),
  m_ncen(0)
{

  // Prepare for computing ordinary moments
  for (auto& product : m_statistics) {
    if (ordinary(product)) {

      m_instOrd.push_back(std::vector<const real*>());
      m_plotOrdinary.push_back(false);
      m_nameOrdinary.push_back(std::string());
      m_ordFieldName.push_back(ctr::FieldName());

      for (auto& term : product) {
        // Put in starting address of instantaneous variable
        m_instOrd[m_nord].push_back(m_particles +
                                    base.control.termOffset(term.quantity) +
                                    term.field);
        if (term.plot) m_plotOrdinary.back() = true;
        // Put in term name+field
        m_nameOrdinary.back() += m_ordFieldName.back()
                               = ctr::FieldName(term.name, term.field);
      }

      ++m_nord;
    }
  }

  if (m_nord) {
    // Storage for all the required ordinary moments
    // +1 for each thread's 0 as center for ordinary moments
    m_ordinary = base.memory.newEntry<real>(m_nthread*(m_nord+1),
                                            REAL,
                                            SCALAR,
                                            "ordinary moments");

    // Put in zero as index of center for ordinary moments in central products
    m_ordinary[m_nord] = 0.0;

    // Prepare for computing central moments
    for (auto& product : m_statistics) {
      if (!ordinary(product)) {

        m_instCen.push_back(std::vector<const real*>());
        m_center.push_back(std::vector<const real*>());
        m_nameCentral.push_back(std::string());

        for (auto& term : product) {
          // Put in starting address of instantaneous variable
          m_instCen[m_ncen].push_back(m_particles +
                                      base.control.termOffset(term.quantity) +
                                      term.field);
          // Put in index of center for central, m_nord for ordinary moment
          m_center[m_ncen].push_back(
            m_ordinary + (!isupper(term.name) ? mean(term) : m_nord));
          m_nameCentral.back() += ctr::FieldName(term.name, term.field);
        }

        ++m_ncen;
      }
    }

    if (m_ncen) {
      // Storage for all the required central moments
      m_central = base.memory.newEntry<real>(m_nthread*m_ncen,
                                             REAL,
                                             SCALAR,
                                             "central moments");
    } // if (m_ncen)
  } // if (m_nord)

} // Roll back changes and rethrow on error
  catch (Exception& e) {
    // No need to clean up if exception thrown from base constructor
    if (e.func() == __PRETTY_FUNCTION__) finalize();
    throw;
  }
  catch (std::exception&) {
    finalize();
    throw;
  }
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }

Statistics::~Statistics() noexcept
//******************************************************************************
//  Destructor
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Statistics::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Single exit point, called implicitly from destructor or explicitly
//!          from anywhere else. Exception safety: no-throw guarantee: never
//!          throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_ncen) m_base.memory.freeEntry(m_central);
  if (m_nord) m_base.memory.freeEntry(m_ordinary);
}

bool
Statistics::ordinary(const std::vector<ctr::Term>& product) const
//******************************************************************************
//  Find out whether product only contains ordinary moment terms
//! \param[in]  product   Vector of terms
//! \author J. Bakosi
//******************************************************************************
{
  // If and only if all terms are ordinary, the product is ordinary
  bool ordinary = true;
  for (auto& term : product) {
    if (term.moment == ctr::Moment::CENTRAL)
      ordinary = false;
  }
  return ordinary;
}

int
Statistics::mean(const ctr::Term& term) const
//******************************************************************************
//  Return mean for fluctuation
//! \param[in]  term      Term (a fluctuation) whose mean to search for
//! \author J. Bakosi
//******************************************************************************
{
  int size = m_ordFieldName.size();
  for (int i=0; i<size; ++i) {
    if (m_ordFieldName[i].name == toupper(term.name) &&
        m_ordFieldName[i].field == term.field) {
       return i;
    }
  }

  Throw(ExceptType::FATAL, std::string("Cannot find mean for variable ")+term);
}

bool
Statistics::plotOrdinary(const int m) const
//******************************************************************************
//  Find out whether ordinary moment is to be plotted
//! \param[in]  m         Moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(m < m_nord, ExceptType::FATAL,
         "Request for an unavailable ordinary moment");
  return m_plotOrdinary[m];
}

const std::string&
Statistics::nameOrdinary(const int m) const
//******************************************************************************
//  Return the name of ordinary moment
//! \param[in]  m         Ordinary-moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(m < m_nord, ExceptType::FATAL,
         "Request for an unavailable ordinary moment");
  return m_nameOrdinary[m];
}

const std::string&
Statistics::nameCentral(const int m) const
//******************************************************************************
//  Return the name of central moment
//! \param[in]  m         Central-moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(m < m_ncen, ExceptType::FATAL,
         "Request for an unavailable central moment");
  return m_nameCentral[m];
}

void
Statistics::estimateOrdinary()
//******************************************************************************
//  Estimate ordinary moments
//! \author J. Bakosi
//******************************************************************************
{
  uint64_t p;
  int tid, i;
  size_t s, j;
  real prod;

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p, i, s, j, prod)
  #endif
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    // Zero ordinary moment accumulators
    memset(m_ordinary + tid*(m_nord+1), 0, m_nord*sizeof(real));

    // Accumulate ordinary moments
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      for (i=0; i<m_nord; ++i) {
        prod = *(m_instOrd[i][0] + p*m_nprop);
        s = m_instOrd[i].size();
        for (j=1; j<s; ++j) prod *= *(m_instOrd[i][j] + p*m_nprop);
        m_ordinary[tid*(m_nord+1) + i] += prod;
      }
    }
  } // omp parallel

  // Collect ordinary moments from all threads
  for (p=1; p<m_nthread; ++p) {
    for (i=0; i<m_nord; ++i) {
      m_ordinary[i] += m_ordinary[p*(m_nord+1) + i];
    }
  }

  // Finish computing ordinary moments
  for (i=0; i<m_nord; ++i) {
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
  uint64_t p;
  int tid, i;
  size_t s, j;
  real prod;

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p, i, s, j, prod)
  #endif
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    // Zero central moment accumulators
    memset(m_central + tid*m_ncen, 0, m_ncen*sizeof(real));

    // Accumulate central moments
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      for (i=0; i<m_ncen; ++i) {
        prod = *(m_instCen[i][0] + p*m_nprop);
        s = m_instCen[i].size();
        for (j=1; j<s; ++j) {
          prod *= *(m_instCen[i][j] + p*m_nprop) - *(m_center[i][j]);
        }
        m_central[tid*m_ncen + i] += prod;
      }
    }
  } // omp parallel

  // Collect central moments from all threads
  for (p=1; p<m_nthread; ++p) {
    for (i=0; i<m_ncen; ++i) {
      m_central[i] += m_central[p*m_ncen + i];
    }
  }

  // Finish computing central moments
  for (i=0; i<m_ncen; ++i) {
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
