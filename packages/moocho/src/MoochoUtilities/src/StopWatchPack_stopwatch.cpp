// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

//#include <iostream>
//#include <iosfwd>

#include "StopWatchPack_stopwatch.hpp"
#include "Teuchos_Time.hpp"

double StopWatchPack::seconds(void)
{
  return Teuchos::Time::wallTime();
}

/*

#ifndef _INTEL_CXX

// Implementation using C standard library.
// In MS VC++ 6.0 the precision is only about 0.05 sec.

#include <time.h>

double StopWatchPack::seconds(void)
{
    static const double secs_per_tick = ((double)1.0) / CLOCKS_PER_SEC;
  const clock_t ticks = clock();
  const double sec = ( (double) ticks ) * secs_per_tick;
  //std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#else	// _INTEL_CXX implementation

// Windows implementation.

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <assert.h>

namespace {

bool seconds_initialized = false;
LARGE_INTEGER start_count, count_freq;	// counts per sec.

inline void seconds_initialize() {
  if( seconds_initialized ) return;
  // Figure out how often the performance counter increments
  ::QueryPerformanceFrequency( &count_freq );
  // Set this thread's priority as high as reasonably possible to prevent
    // timeslice interruptions
    ::SetThreadPriority( ::GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
  // Get the first count.
  TEUCHOS_TEST_FOR_EXCEPT( !(  QueryPerformanceCounter( &start_count )  ) );
  seconds_initialized = true;
}

}	// end namespace

double StopWatchPack::seconds(void)
{
  seconds_initialize();
  LARGE_INTEGER count;
  QueryPerformanceCounter( &count );
  // "QuadPart" is a 64 bit integer (__int64).  VC++ supports them!
  const double
    sec = (double)( count.QuadPart - start_count.QuadPart ) / count_freq.QuadPart;
  //std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
}

#endif	// _INTEL_CXX

*/
