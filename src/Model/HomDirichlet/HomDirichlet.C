//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.C
  \author    J. Bakosi
  \date      Fri Nov 16 09:07:34 2012
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************

#include <iostream>
#include <limits>
#include <cstring>
#include <cmath>

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
  Model(memory, paradigm, "Homogeneous Dirichlet"), m_nscalar(nscalar),
  m_npar(npar)
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
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);

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
  // Initialize the Dirichlet mix model with N-peak delta
  initNpeakDelta();
}

void
HomDirichlet::initNpeakDelta()
//******************************************************************************
//  Initialize with N-peak delta
//! \author  J. Bakosi
//******************************************************************************
{
  // Get MKL VSL stream state pointers
  const VSLStreamStatePtr* str = m_random->getStr(m_rndStr);
  // Get pointer to scalars
  real* scalar = m_memory->getPtr<real>(m_scalar);
  // Precompute number of non-constrained scalars
  int n = m_nscalar-1;

  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    bool accept = false;
    while (!accept) {
      // Generate non-constrained scalars
      real r[n];
      m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD, str[0], n, r, 0.0, 1.0);

      // Compute their sum
      real sum = r[0];
      for (int i=1; i<n; ++i) sum += r[i];

      // Accept if sum is less then 1.0
      if (sum < 1.0) {
        int pN = p*m_nscalar;
        memcpy(scalar+pN, r, n*sizeof(real));   // put in non-constrained ones
        scalar[pN+n] = 1.0-sum;         // the last one is 1.0-(sum of the rest)
        accept = true;
      }
    }

  }

  // Check if sample space is valid
  for (int p=0; p<m_npar; ++p) {
    int pN = p*m_nscalar;
    real sum = scalar[pN];
    for (int i=1; i<m_nscalar; ++i) sum += scalar[pN+i];
    if (fabs(sum-1.0) > numeric_limits<real>::epsilon()) {
      cout << "!";
    }
  }
}
