//******************************************************************************
/*!
  \file      src/Physics/SimplifiedLangevin/SimplifiedLangevin.C
  \author    J. Bakosi
  \date      Sun 20 Jan 2013 01:15:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************

#include <iostream>

#include <sys/time.h>
#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <SimplifiedLangevin.h>

using namespace Quinoa;

SimplifiedLangevin::SimplifiedLangevin(Memory* memory,
                           Paradigm* paradigm,
                           const int& npar,
                           const real time,
                           const int echo,
                           const int nstep) :
  Physics(memory, paradigm, "Simplified Langevin", time, echo, nstep),
  m_npar(npar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  npar     Number of particles
//! \param[in]  time     Maximum time to simulate
//! \param[in]  echo     One-line info in every few time step
//! \param[in]  nstep    Maximum number of time steps to take
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);
}

SimplifiedLangevin::~SimplifiedLangevin()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    //m_memory->freeEntry(m_MEscalar);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in SimplifiedLangevin destructor" << endl; }
#endif // NDEBUG

  if (m_random) { delete m_random; m_random = nullptr; }
}

void
SimplifiedLangevin::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SimplifiedLangevin::echo()
//******************************************************************************
//  Echo informaion on the physics
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Physics: " << m_name << endl;

  cout << endl;
}

void
SimplifiedLangevin::init()
//******************************************************************************
//  Initialize the physics
//! \author  J. Bakosi
//******************************************************************************
{
}
