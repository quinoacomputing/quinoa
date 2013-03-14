//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Wed 13 Mar 2013 08:29:55 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************

#include <cstring>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Statistics.h>
#include <StatException.h>
#include <Paradigm.h>
#include <Control.h>
#include <Mix.h>

using namespace Quinoa;

Statistics::Statistics(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control,
                       Mix* const mix) :
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_nthread(paradigm->nthread()),
  m_npar(control->get<control::NPAR>()),
  m_mix(mix),
  m_nscalar(mix->nscalar()),
  m_statistics(control->get<control::STATISTICS>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object
//! \param[in]  paradigm Parallel programming object
//! \param[in]  control  Control object
//! \param[in]  mix      Mix model object
//! \author  J. Bakosi
//******************************************************************************
{
  // Setup ordinary moments
  m_nord = 0;
  for (auto& product : m_statistics) {
    if (isOrdinary(product)) {
      m_instantaneous.push_back(vector<const real*>());
      m_plotOrdinary.push_back(false);
      m_nameOrdinary.push_back("");
      for (auto& term : product) {
        m_instantaneous[m_nord].push_back(m_mix->scalars() + term.field);
        if (term.plot) {
          m_plotOrdinary.back() = true;
          m_nameOrdinary.back() = term.name;
        }
      }
      ++m_nord;
    }
  }

  // Storage fo all the required ordinary scalar moments
  m_ordinary =
    m_memory->newEntry<real>(m_nthread*m_nord,
                             REAL,
                             SCALAR,
                             "ordinary moments");
}

Statistics::~Statistics()
//******************************************************************************
//  Destructor
//! \author J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    m_memory->freeEntry(m_ordinary);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in Statistics destructor" << endl; }
#endif // NDEBUG
}

bool
Statistics::isOrdinary(const vector<control::Term>& product)
//******************************************************************************
//  Find out whether product only contains ordinary moment terms
//! \author J. Bakosi
//******************************************************************************
{
  // If and only if all terms are ordinary, the product is ordinary
  bool ordinary = true;
  for (auto& term : product) {
    if (term.moment == control::CENTRAL)
      ordinary = false;
  }
  return ordinary;
}

bool
Statistics::plotOrdinary(const int m) const
//******************************************************************************
//  Find out whether ordinary moment is to be plotted
//! \param[in]  m         Moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(m < m_nord, StatException,FATAL,STATEXCEPT_NO_SUCH_MOMENT);
  return m_plotOrdinary[m];
}

string
Statistics::nameOrdinary(const int m) const
//******************************************************************************
//  Return the name of ordinary moment
//! \param[in]  m         Moment index
//! \author J. Bakosi
//******************************************************************************
{
  Assert(m < m_nord, StatException,FATAL,STATEXCEPT_NO_SUCH_MOMENT);
  return m_nameOrdinary[m];
}

void
Statistics::accumulate()
//******************************************************************************
//  Acumulate statistics
//! \author J. Bakosi
//******************************************************************************
{
  int tid, p, i;
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

    // zero ordinary moment accumulators
    memset(m_ordinary + tid*m_nord, 0, m_nord*sizeof(real));

    // accumulate ordinary moments
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      for (i=0; i<m_nord; ++i) {
        prod = *(m_instantaneous[i][0] + p*m_nscalar);
        s = m_instantaneous[i].size();
        for (j=1; j<s; ++j) prod *= *(m_instantaneous[i][j] + p*m_nscalar);
        m_ordinary[tid*m_nord + i] += prod;
      }
    }
  } // omp parallel

  // collect ordinary moments from all threads
  for (p=1; p<m_nthread; ++p) {
    for (i=0; i<m_nord; ++i) {
      m_ordinary[i] += m_ordinary[p*m_nord + i];
    }
  }

  // finish computing ordinary moments
  for (i=0; i<m_nord; ++i) {
    m_ordinary[i] /= m_npar;
  }
}
