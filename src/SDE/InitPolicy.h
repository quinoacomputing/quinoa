//******************************************************************************
/*!
  \file      src/SDE/InitPolicy.h
  \author    J. Bakosi
  \date      Mon 27 Jan 2014 03:31:56 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Initialization policies
  \details   Initialization policies
*/
//******************************************************************************
#ifndef InitPolicy_h
#define InitPolicy_h

#include <cstring>

#include <Types.h>

namespace quinoa {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {
  InitRaw( std::string& policy, const ParProps& particles, uint64_t npar,
           int nprop, int offset, int ncomp, int nthreads )
  {
    policy = "raw";
  }
};

//! Zero initialization policy: zero particle properties
struct InitZero {
  InitZero( std::string& policy, const ParProps& particles, uint64_t npar,
            int nprop, int offset, int ncomp, int nthreads )
  {
    policy = "zero";

    tk::real* ptr = particles.ptr();
    uint64_t size = particles.size();

    // Compute chunk size
    uint64_t chunk = size/nthreads;

    // Each processor zeros its own portion
    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
      #ifdef _OPENMP
      int myid = omp_get_thread_num();
      #else
      int myid = 0;
      #endif

      memset( ptr + chunk*myid, 0, chunk*sizeof(tk::real) );
    }

    // Zero remaining portion
    memset( ptr + chunk*nthreads, 0, (size % nthreads)*sizeof(tk::real) );
  }
};

} // quinoa::

#endif // InitPolicy_h
