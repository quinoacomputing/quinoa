//******************************************************************************
/*!
  \file      src/Base/MKLRandom.C
  \author    J. Bakosi
  \date      Mon 15 Oct 2012 10:03:32 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <iostream>

#include <mkl_vsl.h>

#include <MKLRandom.h>
#include <MemoryException.h>
#include <MKLException.h>

using namespace Quinoa;

MKLRandom::~MKLRandom()
//******************************************************************************
//  Destructor
//! \details Destroy random number generator streams
//! \author  J. Bakosi
//******************************************************************************
{
  typedef vector<RndStreams>::size_type ST;

  try {
    ST tabs = m_table.size();
    for (ST i=0; i<tabs; ++i) {
      // Get pointer to RndStreams
      const RndStreams* s = &m_table[i];
      // Get pointer to array of thread-stream pointers
      VSLStreamStatePtr* stream = s->stream;
      // Delete all thread streams
      for (Int t=0; t<m_nthreads; ++t) {
        if (stream[t] != nullptr &&
            vslDeleteStream(&stream[t]) != VSL_STATUS_OK) {
          cerr << "WARNING: Failed to delete MKL VSL stream" << endl;
        }
      }
      // Free random number table
      m_memory->freeEntry(s->rnd);
      // Delete all thread-stream pointers
      delete [] stream;
    }
  } catch (...) {
    cerr << "WARNING: Exception in MKLRandom::~MKLRandom" << endl;
  }
  // Delete stream tables
  m_table.clear();
}

void
MKLRandom::newStream(VSLStreamStatePtr* stream,
                     const int& brng,
                     const unsigned int& seed)
//******************************************************************************
//  Call MKL's vslNewStream() and handle error
//! \param[out]  stream  VSL stream state descriptor
//! \param[in]   brng    Index of the basic generator to initialize the stream  
//! \param[in]   seed    Initial condition of the stream  
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vslNewStream(stream, brng, seed);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::copyStream(VSLStreamStatePtr* newstream,
                      const VSLStreamStatePtr& srcstream)
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[out]  newstream  Copied random stream descriptor
//! \param[in]   srcstream  Pointer to the stream state structure to be copied
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vslCopyStream(newstream, srcstream);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::skipAheadStream(VSLStreamStatePtr& stream,
                           const long long int& nskip)
//******************************************************************************
//  Call MKL's vslCopyStream() and handle error
//! \param[in]  stream  Pointer to the stream state structure to which the
//!                     block-splitting method is applied
//! \param[in]  nskip   Number of skipped elements
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vslSkipAheadStream(stream, nskip);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::uniform(const Int& method,
                   VSLStreamStatePtr& stream,
                   const Int& n,
                   Real* r,
                   const Real& a,
                   const Real& b)
//******************************************************************************
//  Call MKL's vdRngUniform() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] rnd      Vector of n uniform random numbers
//! \param[in]  a        Left bound
//! \param[in]  b        Right bound
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vdRngUniform(method, stream, n, r, a, b);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::gaussian(const Int& method,
                    VSLStreamStatePtr& stream,
                    const Int& n,
                    Real* r,
                    const Real& a,
                    const Real& b)
//******************************************************************************
//  Call MKL's vdRngGaussian() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] rnd      Vector of n Gaussian random numbers
//! \param[in]  a        mean
//! \param[in]  b        standard deviation
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vdRngGaussian(method, stream, n, r, a, b);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::gamma(const Int& method,
                 VSLStreamStatePtr& stream,
                 const Int& n,
                 Real* r,
                 const Real& alpha,
                 const Real& a,
                 const Real& beta)
//******************************************************************************
//  Call MKL's vdRngGamma() and handle error
//! \param[in]  method   Generation method
//! \param[in]  stream   Pointer to the stream state structure
//! \param[in]  n        Number of random values to be generated
//! \param[out] rnd      Vector of n gamma-distributed random numbers
//! \param[in]  alpha    Shape
//! \param[in]  a        Displacement
//! \param[in]  beta     Scale factor
//! \author  J. Bakosi
//******************************************************************************
{
  Int vslerr = vdRngGamma(method, stream, n, r, alpha, a, beta);
  if (vslerr != VSL_STATUS_OK) throw MKLException(FATAL, vslerr);
}

void
MKLRandom::addTable(const int brng,
                    const Distribution dist,
                    const int method,
                    const long long int number,
                    const string name)
//******************************************************************************
//  Add random number table
//! \param[in]  brng     Basic VSL generator type
//! \param[in]  dist     Type of distribution
//! \param[in]  method   Generation method
//! \param[in]  number   Number of random numbers in table
//! \param[in]  name     Name of random number table
//! \details Create a new random number table. These tables are used to generate
//!          a fix number of fixed property random numbers into tables by
//!          several threads using block-splitting from independent streams
//!          created by VSL's vslskipAheadStream() function.
//! \author  J. Bakosi
//******************************************************************************
{
  // Allocate memory for array of streams of random numbers for several threads
  try {
    m_table.push_back(
      RndStreams(dist,                                 // distribution
                 method,                               // generation method
                 number / m_nthreads,                  // chunk per thread
                 number % m_nthreads,                  // remainder
                 new VSLStreamStatePtr [m_nthreads](), // initialize to zero
                 m_memory->newEntry(number, REAL, SCALAR, name)));
  } catch (bad_alloc&) { throw MemoryException(FATAL, BAD_ALLOC); }

  // Get pointer to newly created RndStreams
  const RndStreams* s = &m_table.back();
  // Get its chunk size
  const long long int chunk = s->chunk;
  // Get pointer to newly created array of thread-stream pointers
  VSLStreamStatePtr* stream = s->stream;

  // Initialize first thread-stream for given distribution using seed
  newStream(&stream[0], brng, m_seed);
  // Initialize the rest of the thread-streams for block-splitting
  for (Int t=0; t<m_nthreads-1; t++) {
    copyStream(&stream[t+1], stream[t]);
    skipAheadStream(stream[t+1], chunk);
  }
}

void
MKLRandom::regenTable(const Int table)
//******************************************************************************
//  Regenerate random numbers in a table
//! \param[in]  table     Table to regenerate
//! \author  J. Bakosi
//******************************************************************************
{
  // Unpack
  const RndStreams* s = &m_table[table];
  const Distribution dist = s->dist;
  const int method = s->method;
  const long long int chunk = s->chunk;
  const long long int remainder = s->remainder;
  VSLStreamStatePtr* stream = s->stream;
  Real* rnd = m_memory->getPtr<Real>(s->rnd);

  // Regenerate
  switch (dist) {
    case UNIFORM:
      #ifdef _OPENMP
      //#pragma omp parallel for
      #endif
      for (Int t=0; t<m_nthreads; ++t) {
        uniform(method, stream[t], chunk, rnd + t*chunk,
                UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);
      }
      uniform(method, stream[0], remainder, rnd + m_nthreads*chunk,
              UNIFORM_LEFT_BOUND, UNIFORM_RIGHT_BOUND);
      break;
    case GAUSSIAN:
      #ifdef _OPENMP
      //#pragma omp parallel for
      #endif
      for (Int t=0; t<m_nthreads; ++t) {
        gaussian(method, stream[t], chunk, rnd + t*chunk,
                 GAUSSIAN_MEAN, GAUSSIAN_STD);
      }
      gaussian(method, stream[0], remainder, rnd + m_nthreads*chunk,
               GAUSSIAN_MEAN, GAUSSIAN_STD);
      break;
    case GAMMA:
      #ifdef _OPENMP
      //#pragma omp parallel for
      #endif
      for (Int t=0; t<m_nthreads; ++t) {
        gamma(method, stream[t], chunk, rnd + t*chunk,
              GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);
      }
      gamma(method, stream[0], remainder, rnd + m_nthreads*chunk,
            GAMMA_SHAPE, GAMMA_DISPLACEMENT, GAMMA_SCALE);
      break;
    default:
      throw MKLException(WARNING, MKLEXCEPT_UNKNOWN_METHOD);
  }
}

void
MKLRandom::regenAllTables()
//******************************************************************************
//  Regenerate random numbers in all tables
//! \author  J. Bakosi
//******************************************************************************
{
  typedef vector<RndStreams>::size_type ST;

  ST tabs = m_table.size();
  for (ST i=0; i<tabs; ++i)
    regenTable(i);
}

// void preprng_streams(int nthreads, VSLStreamStatePtr **stream
//                      #ifndef WALLFUNCTIONS
//                      , int restarted, int samenthreads
//                      #endif
//                      )
// // -----------------------------------------------------------------------------
// // Routine: preprng_streams - Initialize random number generator streams
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // These streams are used to sample a few random numbers at a time
// // with no restrictions on the distribution parameters; prepared for parallel
// // execution.
// // -----------------------------------------------------------------------------
// {
//   #ifndef WALLFUNCTIONS
//   char filename[STRLEN];
//   #endif
// 
//   // allocate memory for array of streams for 'nthreads' threads
//   if(!(*stream=(VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))))
//     ERR("Can't allocate memory!");
// 
//   // initialize stream using the leapfrog technique
//   for (int k=0; k<nthreads; k++) {
//     #ifndef WALLFUNCTIONS
//     if ( restarted && samenthreads )
//     {
//       // construct filename for stream used for a few numbers at a time
//       sprintf( filename, "%s.f.%d", RESTART_FILENAME, k );
//       // load stream
//       CheckVslError( vslLoadStreamF(&(*stream)[k], filename) );
//     }
//     else
//       CheckVslError( vslNewStream(&(*stream)[k], BRNG_FEW, SEED) );
//     #else
//       CheckVslError( vslNewStream(&(*stream)[k], BRNG_FEW, SEED) );
//     #endif
// 
//     CheckVslError( vslLeapfrogStream((*stream)[k], k, nthreads) );
//   }
// }
// 
// void destroyrng_streams( int nthreads, VSLStreamStatePtr **stream )
// // -----------------------------------------------------------------------------
// // Routine: destroyrng_streams - Destroy random number streams
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   // destroy streams
//   for (int k=0; k<nthreads; k++)
//     CheckVslError( vslDeleteStream(&(*stream)[k]) );
// 
//   // pointer to stream
//   free( *stream );
// }
// 
// void saverng_streams( int nthreads
//                       #ifndef WALLFUNCTIONS
//                       , VSLStreamStatePtr *stream
//                       #endif
//                      )
// // -----------------------------------------------------------------------------
// // Routine: saverng_streams - Save the state of random number streams into files
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   char filename[STRLEN];
// 
//   // save the state of random number streams into files
//   for (int k=0; k<nthreads; k++) {
//     #ifndef WALLFUNCTIONS
//     // construct filename for uniform stream used in tables
//     sprintf( filename, "%s.u.%d", RESTART_FILENAME, k );
//     // save stream
//     CheckVslError( vslSaveStreamF(_ustream[k], filename) );
//     #endif
//     
//     // construct filename for Gaussian stream used in tables
//     sprintf( filename, "%s.g.%d", RESTART_FILENAME, k );
//     // save stream
//     CheckVslError( vslSaveStreamF(_gstream[k], filename) );
//     
//     #ifndef WALLFUNCTIONS
//     // construct filename for stream used for a few numbers at a time
//     sprintf( filename, "%s.f.%d", RESTART_FILENAME, k );
//     // save stream
//     CheckVslError( vslSaveStreamF(stream[k], filename) );
//     #endif
//   }
// }
