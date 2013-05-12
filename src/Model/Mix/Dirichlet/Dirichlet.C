//******************************************************************************
/*!
  \file      src/Model/Mix/Dirichlet/Dirichlet.C
  \author    J. Bakosi
  \date      Sun 12 May 2013 06:56:34 AM MDT
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

using namespace std;
using namespace Quinoa;

Dirichlet::Dirichlet(Memory* const memory,
                     Paradigm* const paradigm,
                     Control* const control,
                     real* const scalars)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \param[in]  scalars  Pointer to particle scalars
//! \author  J. Bakosi
//******************************************************************************
try :
  Mix<Dirichlet>(memory,
                 paradigm,
                 control,
                 control->get<control::NSCALAR>(),
                 scalars),
  m_b(control->get<control::B>()),
  m_S(control->get<control::S>()),
  m_k(control->get<control::KAPPA>()),
  m_str(nullptr),
  m_random(nullptr),
  m_rndStr(nullptr)
{

  ErrChk(m_b.size() == static_cast<unsigned int>(m_nscalar), FATAL,
         "Wrong number of Dirichlet model parameters 'b'");
  ErrChk(m_S.size() == static_cast<unsigned int>(m_nscalar), FATAL,
         "Wrong number of Dirichlet model parameters 'S'");
  ErrChk(m_k.size() == static_cast<unsigned int>(m_nscalar), FATAL,
         "Wrong number of Dirichlet model parameters 'k'");

  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(memory, paradigm);
  ErrChk(m_random != nullptr, FATAL,
         "Cannot allocate memory for random number generator");

  // Create random number leapfrog stream
  m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
  // Get array of MKL VSL stream state pointers right away
  m_str = m_random->getStr(m_rndStr);

} // Roll back changes and rethrow on error
  catch (Exception& e) {
    // No need to clean up if exception thrown from base constructor
    if (e.func() == __PRETTY_FUNCTION__) finalize();
    throw;
  }
  catch (exception&) {
    finalize();
    throw;
  }
  catch (...) {
    finalize();
    Throw(UNCAUGHT, "Non-standard exception");
  }

Dirichlet::~Dirichlet() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
Dirichlet::finalize() noexcept
//******************************************************************************
//  Finalize
//! \details Single exit point, called implicitly from destructor or explicitly
//!          from anywhere else. Exception safety: no-throw guarantee: never
//!          throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_random) { delete m_random; m_random = nullptr; }  
}

void
Dirichlet::echo() const
//******************************************************************************
//  Echo information on Dirichlet model
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Dirichlet" << endl;
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
        memcpy(m_scalars + p*m_nscalar, r, m_nscalar*sizeof(real));
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
  real r[m_nscalar];

  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {
    m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                       m_str[0], m_nscalar, r, 0.0, 1.0);
    memcpy(m_scalars + p*m_nscalar, r, m_nscalar*sizeof(real));
  }
}


void
Dirichlet::advance(const real& dt)
//******************************************************************************
//  Advance particles with the Dirichlet model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  int tid, p, i;
  real yn, d;
  real* y;
  real dW[m_nscalar];

  #ifdef _OPENMP
  #pragma omp parallel private(tid, p, i, y, yn, dW, d)
  #endif
  {
    #ifdef _OPENMP
    tid = omp_get_thread_num();
    #else
    tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (p=0; p<m_npar; ++p) {
      // Get access to particle scalars
      y = m_scalars + p*m_nscalar;

      // Compute Nth scalar
      yn = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_nscalar; ++i) yn -= y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                         m_str[tid], m_nscalar, dW, 0.0, 1.0);

      // Advance first m_nscalar (K=N-1) scalars
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
    real* y = m_scalars + p*m_nscalar;
    vector<real> v(y, y+m_nscalar);
    jpdf.insert(v);
  }
}
