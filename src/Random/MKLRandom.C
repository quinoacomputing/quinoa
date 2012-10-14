//******************************************************************************
/*!
  \file      src/Base/MKLRandom.C
  \author    J. Bakosi
  \date      Sat 13 Oct 2012 10:51:46 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     MKL-based random number generator
  \details   MKL-based random number generator
*/
//******************************************************************************

#include <iostream>

#include <mkl_vsl.h>
#include <errcheck.inc>

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
  typedef vector<vector<VSLStreamStatePtr>>::size_type ST;

  ST tabs = table.size();
  for (ST i=0; i<tabs; ++i) {
    // Get pointer to array of thread-stream pointers
    VSLStreamStatePtr* tab = table[i];
    // Delete all thread streams
    for (Int t=0; t<m_nthreads; ++t) {
      if (tab[t] != nullptr && vslDeleteStream(&tab[t]) != VSL_STATUS_OK)
        cerr << "WARNING: Failed to delete MKL VSL stream" << endl;
    }
    // Delete all thread-stream pointers
    delete [] tab;
  }
  // Delete tables
  table.clear();
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
MKLRandom::addTable(Distribution dist, size_t number)
//******************************************************************************
//  Add random number table
//! \author  J. Bakosi
//******************************************************************************
{
  // Allocate memory for array of streams of random numbers for several threads
  try {
    table.push_back(new VSLStreamStatePtr [m_nthreads]()); // initialize to zero
  } catch (bad_alloc&) { throw MemoryException(FATAL, BAD_ALLOC); }

  //double r[10];
  //VSLStreamStatePtr stream;
  //CheckVslError( vslNewStream(&stream, VSL_BRNG_MCG59, m_seed) );
  //CheckVslError( vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, 10,
  //                             r, 0.0, 1.0) );
  //CheckVslError( vslDeleteStream(&stream) );

  // Get pointer to newly created array of stream pointers
  VSLStreamStatePtr* newtab = table.back();

  // Initialize first thread-stream for block-splitting
  newStream(&newtab[0], VSL_BRNG_MCG59, m_seed);

  // Initialize the rest of the thread-streams for block-splitting
  size_t chunk = number / m_nthreads;
  for (Int t=0; t<m_nthreads-1; t++) {
    if (vslCopyStream(&newtab[t+1], newtab[t]) != VSL_STATUS_OK)
      throw MKLException(FATAL, MKL_UNIMPLEMENTED);
    if (vslSkipAheadStream(newtab[t+1], chunk) != VSL_STATUS_OK)
      throw MKLException(FATAL, MKL_UNIMPLEMENTED);
  }
}


// #include <cstdio>
// #include "mkl.h"
// #include "Macros.h"
// #include "Const.h"
// #include "Random.h"
// #include "RandomErrors.h"
// 
// // local data for random number generation in tables
// static int _gchunk, _gremainder;
// static VSLStreamStatePtr *_gstream;
// 
// #ifndef WALLFUNCTIONS
// static int _uchunk, _uremainder;
// static VSLStreamStatePtr *_ustream;
// #endif
// 
// void preprng_tables(int ngr, int nthreads, int restarted, int samenthreads,
//                     #ifndef WALLFUNCTIONS
// 		    int nur, double **ru,
// 		    #endif
// 		    double **rg)
// // -----------------------------------------------------------------------------
// // Routine: preprng_tables - Initialize random number tables
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// // Initialize random number generator streams for parallel generation into
// // tables and allocate memory for random number tables.
// // These streams and tables are used to generate a given (fixed) number of
// // given (fixed property) uniform and Gaussian random numbers into tables in
// // parallel.
// // -----------------------------------------------------------------------------
// {
//   int k;
//   char filename[STRLEN];
// 
// 
//   printf(" * random number table size: %.3g MB\n",
//          (double)(
// 	 #ifndef WALLFUNCTIONS
//          nur+
//          #endif
//          ngr)*sizeof(double)/1024/1024);
//   fflush(stdout);
// 
//   // allocate memory for array of stream of uniform random numbers for
//   // 'nthreads' threads
//   #ifndef WALLFUNCTIONS
//   if(!(_ustream=(VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))))
//     ERR("Can't allocate memory!");
//   #endif
//   if(!(_gstream=(VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))))
//     ERR("Can't allocate memory!");
// 
// 
//   // initialize streams using block-splitting
//   if ( restarted && samenthreads )
//   {
//     #ifndef WALLFUNCTIONS
//     // construct filename for uniform stream used in tables
//     sprintf( filename, "%s.u.0", RESTART_FILENAME );
//     // load stream
//     CheckVslError( vslLoadStreamF(&_ustream[0], filename) );
//     #endif
//     
//     // construct filename for Gaussian stream used in tables
//     sprintf( filename, "%s.g.0", RESTART_FILENAME );
//     // load stream
//     CheckVslError( vslLoadStreamF(&_gstream[0], filename) );
//   }
//   else
//   {
//     #ifndef WALLFUNCTIONS
//     CheckVslError( vslNewStream(&_ustream[0], BRNG_TABLE, SEED) );
//     #endif
//     CheckVslError( vslNewStream(&_gstream[0], BRNG_TABLE, SEED) );
//   }
// 
//   #ifndef WALLFUNCTIONS
//   // compute chunksize and remainder for uniform numbers
//   // ('chunk' random numbers will be generated at once by each processor)
//   _uchunk = nur / nthreads;
//   _uremainder = nur % nthreads;
//   #endif
//   // compute chunksize and remainder for Gaussian numbers
//   // ('chunk' random numbers will be generated at once by each processor)
//   _gchunk = ngr / nthreads;
//   _gremainder = ngr % nthreads;
// 
//   // create SkipAheadStream setting for streams
//   if (restarted && samenthreads)
//     for (k=1; k<nthreads; k++) {
//       #ifndef WALLFUNCTIONS
//       // uniform
//       // construct filename for uniform stream used in tables
//       sprintf( filename, "%s.u.%d", RESTART_FILENAME, k );
//       // load stream
//       CheckVslError( vslLoadStreamF(&_ustream[k], filename) );
//       CheckVslError( vslSkipAheadStream(_ustream[k], _uchunk) );
//       #endif
// 
//       // Gaussian
//       // construct filename for Gaussian stream used in tables
//       sprintf( filename, "%s.g.%d", RESTART_FILENAME, k );
//       // load stream
//       CheckVslError( vslLoadStreamF(&_gstream[k], filename) );
//       CheckVslError( vslSkipAheadStream(_gstream[k], _gchunk) );
//     }
//   else
//     for (k=0; k<nthreads-1; k++) {
//       #ifndef WALLFUNCTIONS
//       // uniform
//       CheckVslError( vslCopyStream(&_ustream[k+1], _ustream[k]) );
//       CheckVslError( vslSkipAheadStream(_ustream[k+1], _uchunk) );
//       #endif
//       // Gaussian
//       CheckVslError( vslCopyStream(&_gstream[k+1], _gstream[k]) );
//       CheckVslError( vslSkipAheadStream(_gstream[k+1], _gchunk) );
//     }
// 
//   #ifndef WALLFUNCTIONS
//   // array to store uniform random numbers
//   if (!(*ru=(double*)malloc(nur*sizeof(double)))) ERR("Can't allocate memory!");
//   #endif
//   // array to store Gaussian random numbers
//   if (!(*rg=(double*)malloc(ngr*sizeof(double)))) ERR("Can't allocate memory!");
// 
//   // initially fill random number tables
//   regenrng_tables(nthreads,
//                   #ifndef WALLFUNCTIONS
//                   *ru,
//                   #endif
//                   *rg );
// }
// 
// void destroyrng_tables(int nthreads,
//                        #ifndef WALLFUNCTIONS
//                        double **ru,
//                        #endif
//                        double **rg )
// // -----------------------------------------------------------------------------
// // Routine: destroyrng_tables - Destroy random number streams and tables
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   // random number tables
//   #ifndef WALLFUNCTIONS
//   free( *ru );
//   #endif
//   free( *rg );
// 
//   // destroy streams
//   for (int k=0; k<nthreads; k++) {
//     #ifndef WALLFUNCTIONS
//     CheckVslError( vslDeleteStream(&_ustream[k]) );
//     #endif
//     CheckVslError( vslDeleteStream(&_gstream[k]) );
//   }
// 
//   // pointers to streams
//   #ifndef WALLFUNCTIONS
//   free( _ustream );
//   #endif
//   free( _gstream );
// }
// 
// void regenrng_tables(int nthreads,
//                      #ifndef WALLFUNCTIONS
//                      double *ru,
//                      #endif
//                      double *rg )
// // -----------------------------------------------------------------------------
// // Routine: regerrng_tables - Regenerate random numbers in tables in parallel
// // Author : J. Bakosi
// // -----------------------------------------------------------------------------
// {
//   // standard uniform between [0 and 1)
//   #ifndef WALLFUNCTIONS
//   #ifdef _OPENMP
//   #pragma omp parallel for
//   #endif
//   for (int k=0; k<nthreads; k++)
//     CheckVslError( vdRngUniform(UNIFORM_METHOD, _ustream[k], _uchunk,
//                                 ru+k*_uchunk, 0.0, 1.0) );
//   // generate remaining portion
//   CheckVslError( vdRngUniform(UNIFORM_METHOD, _ustream[0], _uremainder,
//                               ru+nthreads*_uchunk, 0.0, 1.0) );
//   #endif
// 
//   // Gaussian with zero mean and unit variance
//   #ifdef _OPENMP
//   #pragma omp parallel for
//   #endif
//   for (int k=0; k<nthreads; k++)
//     CheckVslError( vdRngGaussian(GAUSSIAN_METHOD, _gstream[k], _gchunk,
//                                  rg+k*_gchunk, 0.0, 1.0) );
//   // generate remaining portion
//   CheckVslError( vdRngUniform(GAUSSIAN_METHOD, _gstream[0], _gremainder,
//                               rg+nthreads*_gchunk, 0.0, 1.0) );
// }
// 
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
