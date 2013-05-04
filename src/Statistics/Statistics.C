//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Fri 03 May 2013 07:10:19 AM MDT
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

#include <Statistics.h>
#include <Paradigm.h>
#include <Control.h>
#include <Model.h>

using namespace Quinoa;

Statistics::Statistics(Memory* const memory,
                       Paradigm* const paradigm,
                       Control* const control,
                       Model* const model) :
  m_memory(memory),
  m_nthread(paradigm->nthread()),
  m_npar(control->get<control::NPAR>()),
  m_model(model),
  m_nprop(model->nprop()),
  m_statistics(control->get<control::STATISTICS>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object
//! \param[in]  paradigm Parallel programming object
//! \param[in]  control  Control object
//! \param[in]  model    Model objects
//! \author  J. Bakosi
//******************************************************************************
{
  m_nord = 0;
  m_ncen = 0;

  // Prepare for computing ordinary moments
  for (auto& product : m_statistics) {
    if (ordinary(product)) {

      m_instOrd.push_back(vector<const real*>());
      m_plotOrdinary.push_back(false);
      m_nameOrdinary.push_back("");

      for (auto& term : product) {
        // Put in starting address of instantaneous variable
        m_instOrd[m_nord].push_back(m_model->particles() + term.field);
        if (term.plot) m_plotOrdinary.back() = true;
        // Put in term name
        m_nameOrdinary.back() += term.name;
      }

      ++m_nord;
    }
  }

  if (m_nord) {
    // Storage for all the required ordinary moments
    // +1 for each thread's 0 as center for ordinary moments
    m_ordinary = m_memory->newEntry<real>(m_nthread*(m_nord+1),
                                          REAL,
                                          SCALAR,
                                          "ordinary moments");

    // Put in zero as index of center for ordinary moments in central products
    m_ordinary[m_nord] = 0.0;

    // Prepare for computing central moments
    for (auto& product : m_statistics) {
      if (!ordinary(product)) {

        m_instCen.push_back(vector<const real*>());
        m_center.push_back(vector<const real*>());
        m_nameCentral.push_back("");

        for (auto& term : product) {
          // Put in starting address of instantaneous variable
          m_instCen[m_ncen].push_back(m_model->particles() + term.field);
          // Put in index of center for central, m_nord for ordinary moment
          m_center[m_ncen].push_back(
           m_ordinary + (isLower(term.name) ? mean(toUpper(term.name)) : m_nord));
          m_nameCentral.back() += term.name;
        }

        ++m_ncen;
      }
    }

    if (m_ncen) {
      // Storage for all the required central moments
      m_central = m_memory->newEntry<real>(m_nthread*m_ncen,
                                           REAL,
                                           SCALAR,
                                           "central moments");
    } // if (m_ncen)
  } // if (m_nord)
}

Statistics::~Statistics() noexcept
//******************************************************************************
//  Destructor
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author J. Bakosi
//******************************************************************************
{
  // Free memory entries held
  if (m_ncen) m_memory->freeEntry(m_central);
  if (m_nord) m_memory->freeEntry(m_ordinary);
}

bool
Statistics::ordinary(const vector<control::Term>& product)
//******************************************************************************
//  Find out whether product only contains ordinary moment terms
//! \param[in]  product   Vector of terms
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

int
Statistics::mean(const string name) const
//******************************************************************************
//  Return mean for fluctuation
//! \param[in]  name      Name of fluctuation whose mean to search for
//! \author J. Bakosi
//******************************************************************************
{
  int size = m_nameOrdinary.size();
  for (int i=0; i<size; ++i) {
    if (m_nameOrdinary[i] == name) return i;
  }

  throw Exception(FATAL, "Cannot find mean: " + name);
}

string
Statistics::toUpper(const string s) const
//******************************************************************************
//  Convert string to upper case
//! \param[in]  s         String to convert
//! \author J. Bakosi
//******************************************************************************
{
  string upper(s);
  for_each(upper.begin(), upper.end(),
           [](char& c){ c = static_cast<char>(std::toupper(c)); } );
  return upper;
}

bool
Statistics::isLower(const string s) const
//******************************************************************************
//  Return true if string is all lower case
//! \param[in]  s         String to check
//! \author J. Bakosi
//******************************************************************************
{
  bool lower = true;
  for_each(s.begin(), s.end(), [&](char c){ if (isupper(c)) lower = false; } );
  return lower;
}

bool
Statistics::plotOrdinary(const int m) const
//******************************************************************************
//  Find out whether ordinary moment is to be plotted
//! \param[in]  m         Moment index
//! \author J. Bakosi
//******************************************************************************
{
  assert(m < m_nord);
  return m_plotOrdinary[m];
}

const string&
Statistics::nameOrdinary(const int m) const
//******************************************************************************
//  Return the name of ordinary moment
//! \param[in]  m         Ordinary-moment index
//! \author J. Bakosi
//******************************************************************************
{
  assert(m < m_nord);
  return m_nameOrdinary[m];
}

const string&
Statistics::nameCentral(const int m) const
//******************************************************************************
//  Return the name of central moment
//! \param[in]  m         Central-moment index
//! \author J. Bakosi
//******************************************************************************
{
  assert(m < m_ncen);
  return m_nameCentral[m];
}

void
Statistics::estimateOrdinary()
//******************************************************************************
//  Estimate ordinary moments
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
    // Estimate ordinary moments
    estimateOrdinary();

    // Estimate central moments
    if (m_ncen) estimateCentral();
  }
}
