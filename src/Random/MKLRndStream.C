//******************************************************************************
/*!
  \file      src/Random/MKLRndStream.C
  \author    J. Bakosi
  \date      Mon 12 Nov 2012 07:45:28 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation from MKL streams
  \details   Streams are used to generate a few random numbers with no
             restrictions on the distribution parameters using leap-frogging
             between threads
*/
//******************************************************************************

#include <iostream>

#include <MKLRndStream.h>
#include <MKLException.h>

using namespace Quinoa;

MKLRndStream::MKLRndStream(const int nthread,
                           const int brng,
                           const unsigned int seed) :
  m_nthread(nthread)
//******************************************************************************
//  Constructor: Create random number generator leap-frog stream
//! \param[in]  nthread  Number of threads to use
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  seed     Random number generator seed
//! \author  J. Bakosi
//******************************************************************************
{
  Assert(nthread > 0, MKLException,FATAL,MKL_BAD_NTHREADS);

  // Allocate memory for array of stream-pointers for several threads and
  // initialize all to zero
  m_stream = new (nothrow) VSLStreamStatePtr [m_nthread]();
  Assert(m_stream != nullptr, MemoryException,FATAL,BAD_ALLOC);

  // Initialize thread-streams for block-splitting
  for (int t=0; t<m_nthread; ++t) {
    newStream(&m_stream[t], brng, seed);
    leapfrogStream(m_stream[t], t, m_nthread);
  }
}

MKLRndStream::~MKLRndStream()
//******************************************************************************
//  Destructor: Destroy random number generator leap-frog stream
//! \author  J. Bakosi
//******************************************************************************
{
  // Delete all thread streams
  for (int t=0; t<m_nthread; ++t) {
    if (m_stream[t] != nullptr &&
        vslDeleteStream(&m_stream[t]) != VSL_STATUS_OK) {
      cout << "WARNING: Failed to delete MKL VSL stream" << endl;
    }
  }

  // Free all thread-stream pointers
  delete [] m_stream;
}
