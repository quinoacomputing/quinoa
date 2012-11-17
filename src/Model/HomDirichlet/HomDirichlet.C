//******************************************************************************
/*!
  \file      src/Model/HomDirichlet/HomDirichlet.C
  \author    J. Bakosi
  \date      Sat 17 Nov 2012 06:48:12 AM MST
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

#include <Memory.h>
#include <MemoryException.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>
#include <HomDirichlet.h>
#include <Dirichlet.h>

#include <vector>
#include <JPDF.h>
#include <PDFWriter.h>

using namespace Quinoa;

HomDirichlet::HomDirichlet(Memory* memory,
                           Paradigm* paradigm,
                           const int& nscalar,
                           const int& npar,
                           const real time,
                           const int nstep) :
  Model(memory, paradigm, "Homogeneous Dirichlet", time, nstep),
  m_nscalar(nscalar), m_npar(npar)
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  nscalar  Number of mixing scalars
//! \param[in]  npar     Number of particles
//! \param[in]  time     Maximum time to simulate
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
  m_dir = new (nothrow) Dirichlet(nscalar);
  Assert(m_dir != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Allocate memory entry to store the scalars
  m_MEscalar = m_memory->newEntry(npar*nscalar, REAL, SCALAR, "scalar");
  // Get pointer to scalars right away
  m_scalar = m_memory->getPtr<real>(m_MEscalar);
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
    m_memory->freeEntry(m_MEscalar);
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
  initUniform();

//   // Output joint PDF
//   JPDF jpdf(2, 0.1);
//   for (int p=0; p<m_npar; ++p) {
//     int pN = p*m_nscalar;
//     vector<real> v(2,0);
//     v[0] = m_scalar[pN];
//     v[1] = m_scalar[pN+1];
//     jpdf.insert(v);
//   }
//   PDFWriter jpw("../../tmp/jpdf");
//   jpw.write(&jpdf);
}

void
HomDirichlet::initUniform()
//******************************************************************************
//  Initialize scalars with uniform PDF with the last constrained
//! \author  J. Bakosi
//******************************************************************************
{
  // Precompute number of non-constrained scalars
  int n = m_nscalar-1;

  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {

    bool accept = false;
    while (!accept) {
      // Generate non-constrained scalars
      real r[n];
      m_rndStr->uniform(VSL_RNG_METHOD_UNIFORM_STD, m_str[0], n, r, 0.0, 1.0);

      // Compute their sum
      real sum = r[0];
      for (int i=1; i<n; ++i) sum += r[i];

      // Accept if sum is less then 1.0
      if (sum < 1.0) {
        int pN = p*m_nscalar;
        memcpy(m_scalar+pN, r, n*sizeof(real));   // put in non-constrained ones
        m_scalar[pN+n] = 1.0-sum;       // the last one is 1.0-(sum of the rest)
        accept = true;
      }
    }

  }

  // Check if sample space is valid
  for (int p=0; p<m_npar; ++p) {
    int pN = p*m_nscalar;
    real sum = m_scalar[pN];
    for (int i=1; i<m_nscalar; ++i) sum += m_scalar[pN+i];
    if (fabs(sum-1.0) > numeric_limits<real>::epsilon()) {
      cout << "!";
    }
  }
}

void
HomDirichlet::solve()
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
  real dt = 0.1;

  // Time stepping loop
  while (fabs(t-m_time) > numeric_limits<real>::epsilon() && it < m_nstep) {

    // Advance particles
    advance();

    // Echo one-liner info
    report(it, t, dt, hrs2beg, mins2beg, secs2beg, hrs2end, mins2end, secs2end);

    // Increase timestep and iteration counter
    t += dt;
    ++it;
    if (t > m_time) t = m_time;
  }
}

void
HomDirichlet::advance()
//******************************************************************************
//  Advance particles
//! \author  J. Bakosi
//******************************************************************************
{
  for (int p=0; p<m_npar; ++p) {

//    d = kappa[0]*U[pP+1]*(1.0-U[pP+1]-U[pP+2])*dt;
//      if (d > 0.0) d = sqrt(d); else d = 0;//sqrt(-d);
//      U[pP+1] += b[0]/2.0*(S[0]*(1.0-U[pP+1]-U[pP+2]) - (1.0-S[0])*U[pP+1])*dt
//               + d*dW[0];

  }
}

void
HomDirichlet::report(const int it, const real t, const real dt,
                     long int& hrs2beg, long int& mins2beg, long int& secs2beg,
	             long int& hrs2end, long int& mins2end, long int& secs2end)
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
  if (it) secs2end = static_cast<long int>(secs_elapsed*(m_time-t)/(dt*it));
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
