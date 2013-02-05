//******************************************************************************
/*!
  \file      src/Physics/HomogeneousDirichlet/HomogeneousDirichlet.C
  \author    J. Bakosi
  \date      Mon 04 Feb 2013 09:21:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Homogeneous Dirichlet model
  \details   Homogeneous Dirichlet model
*/
//******************************************************************************

#include <iostream>
#include <limits>
#include <cstring>
#include <cmath>

#include <sys/time.h>
#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <HomogeneousDirichlet.h>
#include <Dirichlet.h>

#include <vector>
#include <JPDF.h>
#include <PDFWriter.h>

using namespace Quinoa;

HomogeneousDirichlet::HomogeneousDirichlet(Memory* memory,
                                           Paradigm* paradigm,
                                           const int& nscalar,
                                           const int& npar,
                                           const real time,
                                           const int echo,
                                           const int nstep) :
  Physics(memory, paradigm, "Homogeneous Dirichlet", time, echo, nstep),
  m_N(nscalar), m_npar(npar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
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

  // Instantiate Dirichlet mix model
  m_dir = new (nothrow) Dirichlet(m_N);
  Assert(m_dir != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Allocate memory entry to store the scalars
  m_MEscalar = m_memory->newEntry(npar*m_N, REAL, SCALAR, "scalar");
  // Get pointer to scalars right away
  m_scalar = m_memory->getPtr<real>(m_MEscalar);
}

HomogeneousDirichlet::~HomogeneousDirichlet()
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
    { cout << "WARNING: Exception in HomogeneousDirichlet destructor" << endl; }
#endif // NDEBUG

  if (m_dir) { delete m_dir; m_dir = nullptr; }
  if (m_random) { delete m_random; m_random = nullptr; }
}

void
HomogeneousDirichlet::solve()
//******************************************************************************
//  Solve
//! \author  J. Bakosi
//******************************************************************************
{
  int it=0;
  real t=0.0;
  long int hrs2end=0, mins2end=0, secs2end=0, hrs2beg=0, mins2beg=0, secs2beg=0;

  // Get start time
  gettimeofday(&m_startTime, static_cast<struct timezone*>(0));

  // Set initial time step size
  real dt = 0.05;

  // Time stepping loop
  while (fabs(t-m_time) > numeric_limits<real>::epsilon() && it < m_nstep) {

    // Advance particles
    advance(dt);

    // Echo one-liner info
    if (!(it % m_echo)) {
      report(it, t, dt,
             hrs2beg, mins2beg, secs2beg, hrs2end, mins2end, secs2end);
    }

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_time) t = m_time;
  }

  outJPDF();
}

void
HomogeneousDirichlet::advance(const real dt)
//******************************************************************************
//  Advance particles
//! \author  J. Bakosi
//******************************************************************************
{
  real S[m_N];  S[0] = 5.0/8.0;   S[1] = 2.0/5.0;
  real b[m_N];  b[0] = 0.1;       b[1] = 3.0/2.0;
  real k[m_N];  k[0] = 1.0/80.0;  k[1] = 3.0/10.0;

  int myid, p, i;
  real yn, d;
  real* y;
  real dW[m_N];

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
      y = m_scalar + p*m_N;

      // Compute diagnostic scalar
      yn = 1.0 - y[0];
      #ifdef __INTEL_COMPILER
      #pragma vector always
      #endif
      for (i=1; i<m_N; ++i) yn -= y[i];

      // Generate Gaussian random numbers with zero mean and unit variance
      m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                         m_str[myid], m_N, dW, 0.0, 1.0);

      // Advance prognostic scalars
      for (i=0; i<m_N; ++i) {
        d = k[i]*y[i]*yn*dt;
        if (d > 0.0) d = sqrt(d); else d = 0.0;
        y[i] += b[i]/2.0*(S[i]*yn - (1.0-S[i])*y[i])*dt + d*dW[i];
      }
    } // m_npar
  } // omp parallel
}

void
HomogeneousDirichlet::echo()
//******************************************************************************
//  Echo informaion on homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  cout << "Physics: " << m_name << endl;

  // Echo information on Dirichlet mix model
  m_dir->echo();

  cout << endl;
}

void
HomogeneousDirichlet::init()
//******************************************************************************
//  Initialize homogeneous Dirichlet
//! \author  J. Bakosi
//******************************************************************************
{
  initUniform();
}

void
HomogeneousDirichlet::initUniform()
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
      real r[m_N];
      m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD,
                        m_str[0], m_N, r, 0.0, 1.0);

      // Compute their sum
      real sum = r[0];
      for (int i=1; i<m_N; ++i) sum += r[i];

      // Accept if sum is less then 1.0
      if (sum < 1.0) {
        int pN = p*m_N;
        memcpy(m_scalar+pN, r, m_N*sizeof(real));   // put in scalars
        accept = true;
      }
    }

  }
}

void
HomogeneousDirichlet::initGaussian()
//******************************************************************************
//  Initialize scalars with Gaussian PDF
//! \author  J. Bakosi
//******************************************************************************
{
  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    // Generate scalars
    real r[m_N];
    m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                       m_str[0], m_N, r, 0.0, 1.0);

    int pN = p*m_N;
    memcpy(m_scalar+pN, r, m_N*sizeof(real));   // put in scalars
  }
}

void
HomogeneousDirichlet::outJPDF()
//******************************************************************************
//  Output joint scalar PDF
//! \author  J. Bakosi
//******************************************************************************
{
  JPDF jpdf(2, 0.01);

  for (int p=0; p<m_npar; ++p) {
    real* y = m_scalar + p*m_N;
    vector<real> v(y, y+m_N);
    jpdf.insert(v);
  }

  PDFWriter jpwt(m_memory,"jpdf.txt");
  jpwt.writeTxt(&jpdf);
  PDFWriter jpwg(m_memory,"jpdf.msh");
  jpwg.writeGmsh(&jpdf);
}
