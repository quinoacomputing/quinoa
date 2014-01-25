//******************************************************************************
/*!
  \file      src/SDE/InitPolicy.h
  \author    J. Bakosi
  \date      Fri 24 Jan 2014 07:39:27 AM MST
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
           int nprop, int offset, int ncomp )
  {
    policy = "raw";
  }
};

//! Zero initialization policy: zero particle properties
struct InitZero {
  InitZero( std::string& policy, const ParProps& particles, uint64_t npar,
            int nprop, int offset, int ncomp )
  {
    policy = "zero";
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint64_t p=0; p<npar; ++p) {
      //memset( particles + p*nprop + offset, 0, ncomp*sizeof(tk::real) );
    }
  }
};

} // quinoa::

#endif // InitPolicy_h
