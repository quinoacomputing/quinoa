//  ------------------------------------------------------------------------------------------------------------
//
//  Copyright 2007 Jozsef Bakosi
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  ------------------------------------------------------------------------------------------------------------
//
//  Functions dealing with random numbers. For mor info, see main.cc.
//
//

#include <stdio.h>
#include "mkl.h"
#include "macros.h"
#include "const.h"
#include "random.h"
#include "errcheck.inc"




// local data for random number generation in tables
static int _gchunk, _gremainder;
static VSLStreamStatePtr *_gstream;

#ifndef WALLFUNCTIONS
static int _uchunk, _uremainder;
static VSLStreamStatePtr *_ustream;
#endif





void preprng_tables( int ngr, int nthreads, int restarted, int samenthreads,
                     #ifndef WALLFUNCTIONS
		     int nur, double **ru,
		     #endif
		     double **rg )
//
// initializes random number generator streams for parallel generation into tables,
// allocates memory for random number tables
//
// these streams and tables are used to generate a given (fixed) number of
// given (fixed property) uniform and Gaussian random numbers into tables in parallel
//
{
  int k;
  char filename[STRLEN];


  printf(" * random number table size: %.3g MB\n",
         (double)(
	 #ifndef WALLFUNCTIONS
         nur+
         #endif
         ngr)*sizeof(double)/1024/1024);
  fflush(stdout);

  // allocate memory for array of stream of uniform random numbers for 'nthreads' threads
  #ifndef WALLFUNCTIONS
  if ( !(_ustream = (VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))) )
    ERR("Can't allocate memory!");
  #endif
  if ( !(_gstream = (VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))) )
    ERR("Can't allocate memory!");


  // initialize streams using block-splitting
  if ( restarted && samenthreads )
  {
    #ifndef WALLFUNCTIONS
    // construct filename for uniform stream used in tables
    sprintf( filename, "%s.u.0", RESTART_FILENAME );
    // load stream
    CheckVslError( vslLoadStreamF(&_ustream[0], filename) );
    #endif
    
    // construct filename for Gaussian stream used in tables
    sprintf( filename, "%s.g.0", RESTART_FILENAME );
    // load stream
    CheckVslError( vslLoadStreamF(&_gstream[0], filename) );
  }
  else
  {
    #ifndef WALLFUNCTIONS
    CheckVslError( vslNewStream(&_ustream[0], BRNG_TABLE, SEED) );
    #endif
    CheckVslError( vslNewStream(&_gstream[0], BRNG_TABLE, SEED) );
  }

  #ifndef WALLFUNCTIONS
  // compute chunksize and remainder for uniform numbers
  // ('chunk' random numbers will be generated at once by each processor)
  _uchunk = nur / nthreads;
  _uremainder = nur % nthreads;
  #endif
  // compute chunksize and remainder for Gaussian numbers
  // ('chunk' random numbers will be generated at once by each processor)
  _gchunk = ngr / nthreads;
  _gremainder = ngr % nthreads;

  // create SkipAheadStream setting for streams
  if ( restarted && samenthreads )
    for ( k = 1; k < nthreads; k++ )
    {
      #ifndef WALLFUNCTIONS
      // uniform
      // construct filename for uniform stream used in tables
      sprintf( filename, "%s.u.%d", RESTART_FILENAME, k );
      // load stream
      CheckVslError( vslLoadStreamF(&_ustream[k], filename) );
      CheckVslError( vslSkipAheadStream(_ustream[k], _uchunk) );
      #endif

      // Gaussian
      // construct filename for Gaussian stream used in tables
      sprintf( filename, "%s.g.%d", RESTART_FILENAME, k );
      // load stream
      CheckVslError( vslLoadStreamF(&_gstream[k], filename) );
      CheckVslError( vslSkipAheadStream(_gstream[k], _gchunk) );
    }
  else
    for ( k = 0; k < nthreads-1; k++ )
    {
      #ifndef WALLFUNCTIONS
      // uniform
      CheckVslError( vslCopyStream(&_ustream[k+1], _ustream[k]) );
      CheckVslError( vslSkipAheadStream(_ustream[k+1], _uchunk) );
      #endif
      // Gaussian
      CheckVslError( vslCopyStream(&_gstream[k+1], _gstream[k]) );
      CheckVslError( vslSkipAheadStream(_gstream[k+1], _gchunk) );
    }


  #ifndef WALLFUNCTIONS
  // array to store uniform random numbers
  if ( !(*ru = (double*)malloc(nur*sizeof(double))) ) ERR("Can't allocate memory!");
  #endif
  // array to store Gaussian random numbers
  if ( !(*rg = (double*)malloc(ngr*sizeof(double))) ) ERR("Can't allocate memory!");


  // initially fill random number tables
  regenrng_tables( nthreads,
                   #ifndef WALLFUNCTIONS
		   *ru,
		   #endif
		   *rg );
}







void destroyrng_tables( int nthreads,
                        #ifndef WALLFUNCTIONS
			double **ru,
			#endif
			double **rg )
//
// destroys random number streams and tables
//
{
  int k;


  // random number tables
  #ifndef WALLFUNCTIONS
  free( *ru );
  #endif
  free( *rg );

  // destroy streams
  for ( k = 0; k < nthreads; k++ )
  {
    #ifndef WALLFUNCTIONS
    CheckVslError( vslDeleteStream(&_ustream[k]) );
    #endif
    CheckVslError( vslDeleteStream(&_gstream[k]) );
  }

  // pointers to streams
  #ifndef WALLFUNCTIONS
  free( _ustream );
  #endif
  free( _gstream );

}






void regenrng_tables( int nthreads,
                      #ifndef WALLFUNCTIONS
                      double *ru,
		      #endif
		      double *rg )
//
// regenerates random numbers in tables in parallel
//
{
  int k;


  // standard uniform between [0 and 1)
  #ifndef WALLFUNCTIONS
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( k = 0; k < nthreads; k++ )
    CheckVslError( vdRngUniform(UNIFORM_METHOD, _ustream[k], _uchunk, ru+k*_uchunk, 0.0, 1.0) );
  // generate remaining portion
  CheckVslError( vdRngUniform(UNIFORM_METHOD, _ustream[0], _uremainder, ru+nthreads*_uchunk, 0.0, 1.0) );
  #endif
  
  
  // Gaussian with zero mean and unit variance
  #ifdef _OPENMP
  #pragma omp parallel for
  #endif
  for ( k = 0; k < nthreads; k++ )
    CheckVslError( vdRngGaussian(GAUSSIAN_METHOD, _gstream[k], _gchunk, rg+k*_gchunk, 0.0, 1.0) );
  // generate remaining portion
  CheckVslError( vdRngUniform(GAUSSIAN_METHOD, _gstream[0], _gremainder, rg+nthreads*_gchunk, 0.0, 1.0) );
}








void preprng_streams( int nthreads, VSLStreamStatePtr **stream
                      #ifndef WALLFUNCTIONS
                      , int restarted, int samenthreads
		      #endif
		    )
//
// initializes random number generator streams for parallel generation
//
// this stream is used to sample a few random numbers at a time
// with no restrictions on the distribution parameters
//
// prepared for parallel execution
//
{
  int k;
  #ifndef WALLFUNCTIONS
  char filename[STRLEN];
  #endif


  // allocate memory for array of streams for 'nthreads' threads
  if ( !(*stream = (VSLStreamStatePtr*)malloc(nthreads*sizeof(VSLStreamStatePtr))) )
    ERR("Can't allocate memory!");

  // initialize stream using the leapfrog technique
  for ( k = 0; k < nthreads; k++ )
  {
    #ifndef WALLFUNCTIONS
    if ( restarted && samenthreads )
    {
      // construct filename for stream used for a few numbers at a time
      sprintf( filename, "%s.f.%d", RESTART_FILENAME, k );
      // load stream
      CheckVslError( vslLoadStreamF(&(*stream)[k], filename) );
    }
    else
      CheckVslError( vslNewStream(&(*stream)[k], BRNG_FEW, SEED) );
    #else
      CheckVslError( vslNewStream(&(*stream)[k], BRNG_FEW, SEED) );
    #endif

    CheckVslError( vslLeapfrogStream((*stream)[k], k, nthreads) );
  }
}







void destroyrng_streams( int nthreads, VSLStreamStatePtr **stream )
//
// destroys random number streams
//
{
  int k;


  // destroy streams
  for ( k = 0; k < nthreads; k++ )
    CheckVslError( vslDeleteStream(&(*stream)[k]) );

  // pointer to stream
  free( *stream );
}







void saverng_streams( int nthreads
                      #ifndef WALLFUNCTIONS
                      , VSLStreamStatePtr *stream
                      #endif
                    )
//
// saves the state of all random number streams into files for a later restart
//
{
  int k;
  char filename[STRLEN];


  // save the state of random number streams into files
  for ( k = 0; k < nthreads; k++ )
  {
    #ifndef WALLFUNCTIONS
    // construct filename for uniform stream used in tables
    sprintf( filename, "%s.u.%d", RESTART_FILENAME, k );
    // save stream
    CheckVslError( vslSaveStreamF(_ustream[k], filename) );
    #endif
    
    // construct filename for Gaussian stream used in tables
    sprintf( filename, "%s.g.%d", RESTART_FILENAME, k );
    // save stream
    CheckVslError( vslSaveStreamF(_gstream[k], filename) );
    
    #ifndef WALLFUNCTIONS
    // construct filename for stream used for a few numbers at a time
    sprintf( filename, "%s.f.%d", RESTART_FILENAME, k );
    // save stream
    CheckVslError( vslSaveStreamF(stream[k], filename) );
    #endif
  }
}
