//******************************************************************************
/*!
  \file      src/RNG/MKLRNG.C
  \author    J. Bakosi
  \date      Sat 28 Dec 2013 06:15:57 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <mkl_vsl.h>

#include <MKLRNG.h>
#include <Exception.h>

using quinoa::MKLRNG;

MKLRNG::MKLRNG( int nthreads,
                int brng,
                unsigned int seed,
                int uniform_method,
                int gaussian_method ) :
  m_uniform_method(uniform_method),
  m_gaussian_method(gaussian_method),
  m_nthreads(nthreads)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Throw if not NDEBUG and nthreads invalid
  Assert(nthreads > 0, tk::ExceptType::FATAL, "Need at least one thread");

  // Allocate array of stream-pointers for threads (noexcept)
  m_stream =
    std::unique_ptr< VSLStreamStatePtr[] >( new VSLStreamStatePtr [nthreads] );

  // Initialize thread-streams for block-splitting. These MKL VSL functions
  // dynamically allocate memory, so these calls being in a constructor are a
  // potential memory leak hazard in the presence of exceptions. However,
  // thankfully, the MKL functions below only emit warnings if they encounter
  // errors and always continue. As a result, the constructor finishes, the
  // MKLRNG object gets created, so the destructor will also get called when
  // leaving scope.
  for (auto i=0; i<nthreads; ++i) {
    vslNewStream( &m_stream[i], brng, seed );
    vslLeapfrogStream( m_stream[i], i, nthreads );
  }
}

MKLRNG::~MKLRNG() noexcept
//******************************************************************************
//  Destructor
//! \details No-throw guarantee: this member function never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  // Delete all thread streams
  for (auto i=0; i<m_nthreads; ++i) {
    if (m_stream[i]) {
      vslDeleteStream( &m_stream[i] );
    }
  }
}

void
MKLRNG::uniform(int tid, int num, double* r) const
//******************************************************************************
// Call MKL's uniform RNG
//! \author  J. Bakosi
//******************************************************************************
{
  // method should be exposed to user
  vdRngUniform(m_uniform_method, m_stream[tid], num, r, 0.0, 1.0);
}

void
MKLRNG::gaussian(int tid, int num, double* r) const
//******************************************************************************
// Call MKL's Gaussian RNG
//! \author  J. Bakosi
//******************************************************************************
{
  // method should be exposed to user
  vdRngGaussian(m_gaussian_method, m_stream[tid], num, r, 0.0, 1.0);
}
