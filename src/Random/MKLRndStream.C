//******************************************************************************
/*!
  \file      src/Random/MKLRndStream.C
  \author    J. Bakosi
  \date      Wed May 29 08:34:06 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation from MKL streams
  \details   Streams are used to generate a few random numbers with no
             restrictions on the distribution parameters using leap-frogging
             between threads
*/
//******************************************************************************

#include <iostream>

#include <MKLRndStream.h>
#include <Exception.h>

using namespace Quinoa;

MKLRndStream::MKLRndStream(int nthread,
                           int brng,
                           unsigned int seed)
//******************************************************************************
//  Constructor: Create random number generator leap-frog stream
//! \param[in]  nthread  Number of threads to use
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  seed     Random number generator seed
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the object.
//! \author  J. Bakosi
//******************************************************************************
try :
  m_nthread(nthread),
  m_stream(nullptr)
{

  Assert(nthread > 0, ExceptType::FATAL, "Need at least one thread");

  // Allocate memory for array of stream-pointers for several threads and
  // initialize all to zero
  m_stream = new (nothrow) VSLStreamStatePtr [m_nthread]();
  ErrChk(m_stream != nullptr, ExceptType::FATAL, "Cannot allocate memory");

  // Initialize thread-streams for block-splitting
  for (int t=0; t<m_nthread; ++t) {
    newStream(&m_stream[t], brng, seed);
    leapfrogStream(m_stream[t], t, m_nthread);
  }

} // Roll back changes and rethrow on error
  catch (exception&) {
    finalize();
    throw;
  }
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }


MKLRndStream::~MKLRndStream() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
MKLRndStream::finalize() noexcept
//******************************************************************************
//  Finalize: Destroy random number generator leap-frog stream
//! \details  Single exit point, called implicitly from destructor or explicitly
//!           from anywhere else. Exception safety: no-throw guarantee: never
//!           throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  try {

    // Delete all thread streams
    for (int t=0; t<m_nthread; ++t) {
      if (m_stream[t] != nullptr &&
          vslDeleteStream(&m_stream[t]) != VSL_STATUS_OK) {
        cout << "WARNING: Failed to delete MKL VSL stream" << endl;
      }
    }

  } // emit warning on error
    catch (exception& e) {
      cout << "WARNING: " << e.what() << endl;
    }
    catch (...) {
      cout << "UNKNOWN EXCEPTION in MKLRndStream's destructor" << endl
           << "Continuing anyway..." << endl;
    }

  // Free all thread-stream pointers
  delete [] m_stream;
}
