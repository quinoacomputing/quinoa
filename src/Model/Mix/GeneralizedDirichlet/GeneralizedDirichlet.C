//******************************************************************************
/*!
  \file      src/Model/Mix/GeneralizedDirichlet/GeneralizedDirichlet.C
  \author    J. Bakosi
  \date      Sun 03 Mar 2013 12:47:10 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     The generalized Dirichlet mix model
  \details   The generalized Dirichlet mix model
*/
//******************************************************************************

#include <cstring>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <GeneralizedDirichlet.h>
#include <Mix.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <JPDF.h>
#include <Control.h>
#include <MixException.h>

using namespace std;
using namespace Quinoa;

GeneralizedDirichlet::GeneralizedDirichlet(Memory* const memory,
                                           Paradigm* const paradigm,
                                           Control* const control) :
  Mix(memory, paradigm, control, "Generalized Dirichlet"),
  m_b(control->get<B>()),
  m_S(control->get<S>()),
  m_k(control->get<KAPPA>()),
  m_c(control->get<C>())
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
  Assert(m_c.size() == static_cast<unsigned int>(m_nscalar*(m_nscalar-1)/2),
         MixException, FATAL, BAD_MODEL_PARAMETERS);

  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  Assert(m_random != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

  // Allocate memory entry to store all the scalars
  m_allScalars = m_memory->newEntry<real>(m_npar*m_nscalar,
                                          REAL,
                                          SCALAR,
                                          "allScalars");
}

GeneralizedDirichlet::~GeneralizedDirichlet()
//******************************************************************************
//  Destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    m_memory->freeEntry(m_allScalars);
#ifndef NDEBUG
  } catch (...)
    { cout << "WARNING: Exception in HomMix destructor" << endl; }
#endif // NDEBUG

  if (m_random) { delete m_random; m_random = nullptr; }  
}

void
GeneralizedDirichlet::echo()
//******************************************************************************
//  Echo information on the generalized Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
}

void
GeneralizedDirichlet::init()
//******************************************************************************
//  Initialize scalars
//! \author  J. Bakosi
//******************************************************************************
{
  initUniform();
}

void
GeneralizedDirichlet::initUniform()
//******************************************************************************
//  Initialize scalars with uniform PDF with the last constrained
//! \author  J. Bakosi
//******************************************************************************
{
  real r[m_nscalar];

  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    bool accept = false;
    while (!accept) {
      // Generate scalars
      m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD,
                        m_str[0], m_nscalar, r, 0.0, 1.0);

      // Compute their sum
      real sum = r[0];
      for (int i=1; i<m_nscalar; ++i) sum += r[i];

      // Accept if sum is less then 1.0
      if (sum < 1.0) {
        memcpy(m_allScalars + p*m_nscalar, r, m_nscalar*sizeof(real));
        accept = true;
      }
    }

  }
}

void
GeneralizedDirichlet::advance(const real dt)
//******************************************************************************
//  Advance particles with the generalized Dirichlet model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int myid, p, i, j, k;
  real d, a;
  real* y;
  real dW[m_nscalar];
  real Y[m_nscalar];    //!< Y_i = 1 - sum_{k=1}^{i} y_k
  real U[m_nscalar];    //!< U_i = prod_{j=1}^{m_nscalar-i} 1/sY_{m_nscalar-j}

  #ifdef _OPENMP
  #pragma omp parallel private(myid, p, i, j, k, d, a, y, Y, U, dW)
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
      y = m_allScalars + p*m_nscalar;

      Y[0] = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_nscalar; ++i) Y[i] = Y[i-1] - y[i];

      U[m_nscalar-1] = 1.0;
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=m_nscalar-2; i>=0; --i) U[i] = U[i+1]/Y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                         m_str[myid], m_nscalar, dW, 0.0, 1.0);

      // Advance first m_nscalar (K=N-1) scalars
      k=0;
      a=0.0;
      for (i=0; i<m_nscalar; ++i) {
        d = m_k[i]*y[i]*Y[m_nscalar-1]*U[i]*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        for (j=i; j<m_nscalar-1; ++j) a += m_c[k++]/Y[j];
        y[i] += U[i]/2.0*m_b[i]*
             ((m_S[i]*Y[m_nscalar-1] - (1.0-m_S[i])*y[i]) +
               y[i]*Y[m_nscalar-1]*a)*dt + d*dW[i];
      }
    } // m_npar
  } // omp parallel
}

void
GeneralizedDirichlet::jpdf(JPDF& jpdf)
//******************************************************************************
//  Estimate joint scalar probability density function
//! \author  J. Bakosi
//******************************************************************************
{
  for (int p=0; p<m_npar; ++p) {
    real* y = m_allScalars + p*m_nscalar;
    vector<real> v(y, y+m_nscalar);
    jpdf.insert(v);
  }
}
