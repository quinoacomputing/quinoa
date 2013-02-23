//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Sat 23 Feb 2013 11:57:23 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Dirichlet mix model
  \details   Dirichlet mix model
*/
//******************************************************************************

#include <cstring>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Dirichlet.h>
#include <Mix.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <JPDF.h>
#include <Control.h>
#include <MixException.h>

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control) :
  Mix(memory, paradigm, control, "Dirichlet"),
  m_b(control->get<B>()),
  m_S(control->get<S>()),
  m_k(control->get<KAPPA>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(m_b.size() == static_cast<unsigned int>(m_nscalar),
         MixException, FATAL, BAD_MODEL_PARAMETERS);
  Assert(m_S.size() == static_cast<unsigned int>(m_nscalar),
         MixException, FATAL, BAD_MODEL_PARAMETERS);
  Assert(m_k.size() == static_cast<unsigned int>(m_nscalar),
         MixException, FATAL, BAD_MODEL_PARAMETERS);

  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

  // Allocate memory entry to store the scalars
  m_MEscalar = m_memory->newEntry(m_npar*m_nscalar, REAL, SCALAR, "scalar");
  // Get pointer to scalars right away
  m_scalar = m_memory->getPtr<real>(m_MEscalar);
}

Dirichlet::~Dirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    m_memory->freeEntry(m_MEscalar);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in HomMix destructor" << endl; }
#endif // NDEBUG

  if (m_random) { delete m_random; m_random = nullptr; }  
}

void
Dirichlet::echo()
//******************************************************************************
//  Echo information on Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
}

void
Dirichlet::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
  initUniform();
}

void
Dirichlet::initUniform()
//******************************************************************************
//  Initialize scalars with uniform PDF with the last constrained
//! \author  J. Bakosi
//******************************************************************************
{
  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    bool accept = false;
    while (!accept) {
      // Generate scalars
      real r[m_nscalar];
      m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD,
                        m_str[0], m_nscalar, r, 0.0, 1.0);

      // Compute their sum
      real sum = r[0];
      for (int i=1; i<m_nscalar; ++i) sum += r[i];

      // Accept if sum is less then 1.0
      if (sum < 1.0) {
        int pN = p*m_nscalar;
        memcpy(m_scalar+pN, r, m_nscalar*sizeof(real));   // put in scalars
        accept = true;
      }
    }

  }
}

void
Dirichlet::initGaussian()
//******************************************************************************
//  Initialize scalars with Gaussian PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    // Generate scalars
    real r[m_nscalar];
    m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                       m_str[0], m_nscalar, r, 0.0, 1.0);

    int pN = p*m_nscalar;
    memcpy(m_scalar+pN, r, m_nscalar*sizeof(real));   // put in scalars
  }
}


void
Dirichlet::advance(const real dt)
//******************************************************************************
//  Advance particles with the Dirichlet model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int myid, p, i;
  real yn, d;
  real* y;
  real dW[m_nscalar];

  #ifdef _OPENMP
  #pragma omp parallel private(myid, p, y, yn, i, dW, d)
  #endif // _OPENMP
  {
    #ifdef _OPENMP
    myid = omp_get_thread_num();
    #else
    myid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      // Get access to particle scalars
      y = m_scalar + p*m_nscalar;

      // Compute diagnostic scalar
      yn = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_nscalar; ++i) yn -= y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                         m_str[myid], m_nscalar, dW, 0.0, 1.0);

      // Advance prognostic scalars
      for (i=0; i<m_nscalar; ++i) {
        d = m_k[i]*y[i]*yn*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        y[i] += m_b[i]/2.0*(m_S[i]*yn - (1.0-m_S[i])*y[i])*dt + d*dW[i];
      }
    } // m_npar
  } // omp parallel
}

void
Dirichlet::jpdf(JPDF& jpdf)
//******************************************************************************
//  Estimate joint scalar probability density function
//! \author  J. Bakosi
//******************************************************************************
{
  for (int p=0; p<m_npar; ++p) {
    real* y = m_scalar + p*m_nscalar;
    vector<real> v(y, y+m_nscalar);
    jpdf.insert(v);
  }
}
