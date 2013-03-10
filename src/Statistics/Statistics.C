//******************************************************************************
/*!
  \file      src/Statistics/Statistics.C
  \author    J. Bakosi
  \date      Sun 10 Mar 2013 01:46:34 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Statistics
  \details   Statistics
*/
//******************************************************************************

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
  m_mix(mix),
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
  // Setup means
  for (auto& product : m_statistics) {
    if (product.size() == 1 && product[0].moment == control::ORDINARY) {
      control::Term term = product[0];
      if (term.quantity == control::TRANSPORTED_SCALAR) {
        cout << term.name << ", " << term.field << endl;
        m_instantaneous.push_back(m_mix->scalars() + term.field);
      }
    }
  }
  for (auto& inst : m_instantaneous) cout << inst << endl;

  // Store number of ordinary moments
  m_nord = m_instantaneous.size();

  // Allocate memory to store all the required ordinary scalar moments
  m_ordinary =
    m_memory->newEntry<real>(m_nthread*m_nord,
                             REAL,
                             SCALAR,
                             "ordinary scalar moments");

  // Setup ordinary moment accumulators for each thread
  for (int i=0; i<m_nthread; ++i) {
    m_ordinaryThread.push_back(m_ordinary.ptr + i*m_nord);
  }
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

void
Statistics::accumulate()
//******************************************************************************
//  Acumulate statistics
//! \author J. Bakosi
//******************************************************************************
{
}
