//******************************************************************************
/*!
  \file      src/Model/Hydro/SimplifiedLangevin/SimplifiedLangevin.C
  \author    J. Bakosi
  \date      Mon May  6 13:06:39 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Simplified Langevin hydrodynamics model
  \details   Simplified Langevin hydrodynamics model
*/
//******************************************************************************

#include <iostream>
#include <cstring>
#include <cmath>

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <Macro.h>
#include <SimplifiedLangevin.h>
#include <Hydro.h>
#include <Control.h>
#include <MKLRandom.h>
#include <MKLRndStream.h>

using namespace std;
using namespace Quinoa;

SimplifiedLangevin::SimplifiedLangevin(Memory* const memory,
                                       Paradigm* const paradigm,
                                       Control* const control) :
  Hydro(memory, paradigm, control, "Simplified Langevin"),
  m_C0(control->get<control::C0>())
//******************************************************************************
//  Constructor
//! \param[in]  memory   Memory object pointer
//! \param[in]  paradigm Parallel programming object pointer
//! \param[in]  control  Control object pointer
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate random number generator
  m_random = new (nothrow) MKLRandom(m_memory, m_paradigm);
  if (m_random == nullptr)
    throw Exception(FATAL, "Cannot allocate memory for random number generator "
                           "in SimplifiedLangevin constructor");

  try {

    // Create random number leapfrog stream
    m_rndStr = m_random->addStream(VSL_BRNG_MCG59, 0);
    // Get array of MKL VSL stream state pointers right away
    m_str = m_random->getStr(m_rndStr);

    // Allocate memory to store all the particle properties
    m_particles = m_memory->newEntry<real>(m_npar*m_nprop,
                                           REAL,
                                           SCALAR,
                                           "SLM particles");

  } // roll back changes and rethrow on error
    catch (exception&) {
      finalize();
      throw;
    }
    catch (...) {
      finalize();
      throw Exception(UNCAUGHT);
    }
}

SimplifiedLangevin::~SimplifiedLangevin() noexcept
//******************************************************************************
//  Destructor
//! \details    Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
SimplifiedLangevin::finalize() noexcept
//******************************************************************************
//  Finalize, single exit point called from the destructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Free memory entries held
  m_memory->freeEntry(m_particles);

  if (m_random) { delete m_random; m_random = nullptr; }
}

void
SimplifiedLangevin::echo() const
//******************************************************************************
//  Echo information on the simplified Langevin model
//! \author  J. Bakosi
//******************************************************************************
{
}

void
SimplifiedLangevin::init()
//******************************************************************************
//  Initialize particle properties
//! \author  J. Bakosi
//******************************************************************************
{
  real r[m_nprop];

  // Generate initial values for all scalars for all particles
  for (int p=0; p<m_npar; ++p) {
    m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
                       m_str[0], m_nprop, r, 0.0, 1.0);
    memcpy(m_particles + p*m_nprop, r, m_nprop*sizeof(real));
  }
}

void
SimplifiedLangevin::advance(const real dt)
//******************************************************************************
//  Advance particles with the simplified Langevin model
//! \param[in]  dt   Time step size
//! \author  J. Bakosi
//******************************************************************************
{
  IGNORE(dt);
//   int tid, p, i;
//   real d, tke, eps, S=1.0;
//   real* X;
//   real* U;
//   real dW[3], rs[6];
// 
//   #ifdef _OPENMP
//   //#pragma omp parallel private(tid, p, i, y, yn, dW, d)
//   #endif
//   {
//     #ifdef _OPENMP
//     tid = omp_get_thread_num();
//     #else
//     tid = 0;
//     #endif
// 
//     #ifdef _OPENMP
//     #pragma omp for
//     #endif
//     for (p=0; p<m_npar; ++p) {
//       // Get access to particle position, velocity
//       X = m_particles + p*m_nprop;
//       U = X + 3;
// 
//       // Generate Gaussian random numbers with zero mean and unit variance
//       m_rndStr->gaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER,
//                          m_str[tid], 3, dW, 0.0, 1.0);
// 
//       tke = 0.5*(rs[0] + rs[1] + rs[2]);
//       eps = S*tke/2.36;
// 
//       // Advance velocity
//       d = m_C0*eps*dt;
//       if (d > 0.0) d = sqrt(d); else d = 0.0;
//       //U[0] += ()*dt + d*dW[i];
//     } // m_npar
//   } // omp parallel
}
