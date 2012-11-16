//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.C
  \author    J. Bakosi
  \date      Thu Nov 15 17:36:57 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************

#include <iostream>

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <HomDirichlet.h>
#include <Dirichlet.h>

using namespace Quinoa;

HomDirichlet::HomDirichlet(Memory* memory,
                           Paradigm* paradigm,
                           const int& nscalar,
                           const int& npar) :
  Model(memory, paradigm, "Homogeneous Dirichlet"), m_npar(npar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  npar     Number of particles
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Create random number leapfrog stream
  m_rndStream = m_random->addStream(VSL_BRNG_MCG59, 0);

  // Instantiate Dirichlet mix model
  m_dir = new (nothrow) Dirichlet(nscalar);
  Assert(m_dir != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Allocate memory entry to store the scalars
  m_scalar = m_memory->newEntry(npar*nscalar, REAL, SCALAR, "scalar");
}

HomDirichlet::~HomDirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // No error checking done and no exceptions thrown in debug mode
  try {
#endif // NDEBUG
    m_memory->freeEntry(m_scalar);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in HomDirichlet::~HomDirichlet" << endl; }
#endif // NDEBUG

  if (m_dir) { delete m_dir; m_dir = nullptr; }
  if (m_random) { delete m_random; m_random = nullptr; }
}

void
HomDirichlet::echo()
//******************************************************************************
//  Echo informaion on homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Model: " << m_name << endl;

  // Echo information on Dirichlet mix model
  m_dir->echo();

  cout << endl;
}

void
HomDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  // Initialize the Dirichlet mix model
  m_dir->init(m_npar, m_memory->getPtr<real>(m_scalar));
}
