//******************************************************************************
/*!
  \file      src/Physics/Physics.C
  \author    J. Bakosi
  \date      Wed 20 Feb 2013 09:05:30 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics base
  \details   Physics base
*/
//******************************************************************************

#include <iostream>

#include <sys/time.h>

#include <Memory.h>
#include <Physics.h>
#include <Control.h>

using namespace std;
using namespace Quinoa;
using namespace control;

Physics::Physics(Memory* const memory,
                 Paradigm* const paradigm,
                 Control* const control) :
  m_memory(memory),
  m_paradigm(paradigm),
  m_control(control),
  m_term(control->get<TERM>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \author  J. Bakosi
//******************************************************************************
{
}

void
Physics::report(const int it,
                const real t,
                const real dt,
                long int& hrs2beg,
                long int& mins2beg,
                long int& secs2beg,
	        long int& hrs2end,
                long int& mins2end,
                long int& secs2end)
//******************************************************************************
//  One-liner report
//! \param[in]  it        Iteration counter
//! \param[in]  t         Time counter
//! \param[in]  dt        Time step size
//! \param[in]  hrs2beg   Hours elapsed
//! \param[in]  mins2beg  Minutes elapsed
//! \param[in]  secs2beg  Seconds elapsed
//! \param[in]  hrs2end   Estimated hours until finish
//! \param[in]  mins2end  Estimated minutes until finish
//! \param[in]  secs2end  Estimate seconds until finish
//! \author  J. Bakosi
//******************************************************************************
{
  struct timeval cur_time;
  long int secs_elapsed;

  // Get current time
  gettimeofday( &cur_time, (struct timezone*)0 );
  secs_elapsed = ((cur_time.tv_sec - m_startTime.tv_sec) * 1000000 +
                  (cur_time.tv_usec - m_startTime.tv_usec)) / 1000000;

  // Calculate elapsed time
  secs2beg = secs_elapsed;
  mins2beg = secs2beg/60;
  hrs2beg = mins2beg/60;
  if (secs2beg >= 60) secs2beg %= 60;
  if (secs2beg >= 60) secs2beg %= 60;
  if (mins2beg >= 60) mins2beg %= 60;

  // Estimate time remaining
  if (it) secs2end = static_cast<long int>(secs_elapsed*(m_term-t)/(dt*it));
  else secs2end = 0;
  mins2end = secs2end/60;
  hrs2end = mins2end/60;
  if (secs2end >= 60) secs2end %= 60;
  if (secs2end >= 60) secs2end %= 60;
  if (mins2end >= 60) mins2end %= 60;

  cout << "it = " << it << ", t = " << t << "\t dt = " << dt << "\t"
       << hrs2beg << ":" << mins2beg << ":" << secs2beg << "\t"
       << hrs2end << ":" << mins2end << ":" << secs2end << endl;
}
