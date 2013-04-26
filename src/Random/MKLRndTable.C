//******************************************************************************
/*!
  \file      src/Random/MKLRndTable.C
  \author    J. Bakosi
  \date      Fri Apr 26 16:04:41 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation into tables using Intel's MKL
  \details   Tables are used to generate a fix number of fixed property random
             numbers by several threads using block-splitting.
*/
//******************************************************************************

#include <iostream>

#include <MKLRndTable.h>
#include <MKLException.h>

using namespace Quinoa;

MKLRndTable::MKLRndTable(Memory* const memory,
                         const int nthread,
                         const int brng,
                         const RndDist dist,
                         const int method,
                         const unsigned int seed,
                         const long long int number,
                         const string name) :
  m_memory(memory),
  m_nthread(nthread),
  m_dist(dist),
  m_method(method),
  m_chunk(number / nthread),
  m_remainder(number % nthread)
//******************************************************************************
//  Constructor: Create random number table
//! \param[in]  memory   Memory object pointer
//! \param[in]  nthread  Number of threads to use
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  dist     Type of distribution
//! \param[in]  method   Generation method
//! \param[in]  seed     Random number generator seed
//! \param[in]  number   Total number of random numbers in table
//! \param[in]  name     MemoryEntry name of the random number table
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(nthread > 0, MKLException,FATAL,MKL_BAD_NTHREADS);
  Assert(number > 0, MKLException,FATAL, MKL_BAD_NUMBER);
  Assert(name.size() > 0, MemoryException,FATAL,EMPTY_NAME);

  // Allocate memory for array of stream-pointers for several threads and
  // initialize all to zero
  m_stream = new (nothrow) VSLStreamStatePtr [m_nthread]();
  Assert(m_stream != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Initialize first thread-stream for given distribution using seed
  newStream(&m_stream[0], brng, seed);
  // Initialize the rest of the thread-streams for block-splitting
  for (int t=0; t<m_nthread-1; t++) {
    copyStream(&m_stream[t+1], m_stream[t]);
    skipAheadStream(m_stream[t+1], m_chunk);
  }

  // Allocate array to store random numbers
  m_rnd = m_memory->newEntry<real>(number, REAL, SCALAR, name);
}

MKLRndTable::~MKLRndTable() noexcept
//******************************************************************************
//  Destructor: Destroy random number table
//! \author  J. Bakosi
//******************************************************************************
{
#ifndef NDEBUG  // Error checking and exceptions only in debug mode
  try {
#endif // NDEBUG
    // Delete all thread streams
    for (int t=0; t<m_nthread; ++t) {
      if (m_stream[t] != nullptr &&
          vslDeleteStream(&m_stream[t]) != VSL_STATUS_OK) {
        cout << "WARNING: Failed to delete MKL VSL stream" << endl;
      }
    }
    // Free all thread-stream pointers
    delete [] m_stream;
    // Free array storing random numbers
    m_memory->freeEntry(m_rnd);
#ifndef NDEBUG
  } catch (...) {
    cout << "WARNING: Exception in MKLRndTable::~MKLRndTable" << endl;
  }
#endif // NDEBUG
}


void
MKLRndTable::generate() const
//******************************************************************************
//  Regenerate random numbers in a table
//! \author  J. Bakosi
//******************************************************************************
{
  switch (m_dist) {
    case UNIFORM:
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int t=0; t<m_nthread; ++t) {
        uniform(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
                UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);
      }
      uniform(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
              UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);
      break;
    case GAUSSIAN:
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int t=0; t<m_nthread; ++t) {
        gaussian(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
                 GAUSSIAN_MEAN, GAUSSIAN_STD);
      }
      gaussian(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
               GAUSSIAN_MEAN, GAUSSIAN_STD);
      break;
    case GAMMA:
      #ifdef _OPENMP
      #pragma omp parallel for
      #endif
      for (int t=0; t<m_nthread; ++t) {
        gamma(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
              GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);
      }
      gamma(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
            GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);
      break;
    default:
      Throw(MKLException,WARNING,MKL_UNKNOWN_METHOD);
  }
}
