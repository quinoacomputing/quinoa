//******************************************************************************
/*!
  \file      src/Random/MKLRndTable.C
  \author    J. Bakosi
  \date      Sun 15 Sep 2013 05:52:08 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generation into tables using Intel's MKL
  \details   Tables are used to generate a fix number of fixed property random
             numbers by several threads using block-splitting.
*/
//******************************************************************************

#include <iostream>

#include <MKLRndTable.h>

using namespace quinoa;

MKLRndTable::MKLRndTable(const Memory& memory,
                         int nthread,
                         int brng,
                         RndDist dist,
                         int method,
                         unsigned int seed,
                         long long int number,
                         const std::string& name)
//******************************************************************************
//  Constructor: Create random number table
//! \param[in]  memory   Memory manager
//! \param[in]  nthread  Number of threads to use
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  dist     Type of distribution
//! \param[in]  method   Generation method
//! \param[in]  seed     Random number generator seed
//! \param[in]  number   Total number of random numbers in table
//! \param[in]  name     MemoryEntry name of the random number table
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the object.
//! \author  J. Bakosi
//******************************************************************************
try :
  m_memory(memory),
  m_nthread(nthread),
  m_dist(dist),
  m_method(method),
  m_chunk(number / nthread),
  m_remainder(number % nthread),
  m_stream(nullptr),
  m_rnd()
{

  Assert(nthread > 0, ExceptType::FATAL, "Need at least one thread");
  Assert(number > 0, ExceptType::FATAL,
         "Number of random numbers must be positive");
  Assert(name.size() > 0, ExceptType::FATAL,
         "MKLRndTable's name must be nonempty");

  // Allocate memory for array of stream-pointers for several threads and
  // initialize all to zero
  m_stream = new (std::nothrow) VSLStreamStatePtr [m_nthread]();
  if (m_stream == nullptr)
    Exception(ExceptType::FATAL, "Cannot allocate memory");

  // Initialize first thread-stream for given distribution using seed
  newStream(&m_stream[0], brng, seed);
  // Initialize the rest of the thread-streams for block-splitting
  for (int t=0; t<m_nthread-1; t++) {
    copyStream(&m_stream[t+1], m_stream[t]);
    skipAheadStream(m_stream[t+1], m_chunk);
  }

  // Allocate array to store random numbers
  m_rnd = m_memory.newEntry<real>(number, REAL, SCALAR, name);

} // Roll back changes and rethrow on error
  catch (std::exception&) {
    finalize();
    throw;
  }
  catch (...) {
    finalize();
    Throw(ExceptType::UNCAUGHT, "Non-standard exception");
  }



MKLRndTable::~MKLRndTable() noexcept
//******************************************************************************
//  Destructor
//! \details Exception safety: no-throw guarantee: never throws exceptions.
//! \author  J. Bakosi
//******************************************************************************
{
  finalize();
}

void
MKLRndTable::finalize() noexcept
//******************************************************************************
//  Finalize: Destroy random number table
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
        std::cout << "WARNING: Failed to delete MKL VSL stream" << std::endl;
      }
    }

  } // emit warning on error
    catch (std::exception& e) {
      std::cout << "WARNING: " << e.what() << std::endl;
    }
    catch (...) {
      std::cout << "UNKNOWN EXCEPTION in MKLRndTable's destructor" << std::endl
                << "Continuing anyway..." << std::endl;
    }

  // Free all thread-stream pointers
  delete [] m_stream;
  // Free array storing random numbers
  m_memory.freeEntry(m_rnd);
}


void
MKLRndTable::generate() const
//******************************************************************************
//  Regenerate random numbers in a table
//! \details    Exception safety: strong guarantee: if an exception is thrown,
//!             there are no changes to the object.
//! \author  J. Bakosi
//******************************************************************************
{
  if (m_dist == RndDist::UNIFORM) {

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int t=0; t<m_nthread; ++t) {
      uniform(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
              UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);
    }
    uniform(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
            UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);

  } else if (m_dist == RndDist::GAUSSIAN) {

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int t=0; t<m_nthread; ++t) {
      gaussian(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
               GAUSSIAN_MEAN, GAUSSIAN_STD);
    }
    gaussian(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
             GAUSSIAN_MEAN, GAUSSIAN_STD);

  } else if (m_dist == RndDist::GAMMA) {

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int t=0; t<m_nthread; ++t) {
      gamma(m_method, m_stream[t], m_chunk, m_rnd + t*m_chunk,
            GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);
    }
    gamma(m_method, m_stream[0], m_remainder, m_rnd + m_nthread*m_chunk,
          GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);

  } else {
    Throw(ExceptType::WARNING, "Unknown random number distribution");
  }
}
