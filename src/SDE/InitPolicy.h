//******************************************************************************
/*!
  \file      src/SDE/InitPolicy.h
  \author    J. Bakosi
  \date      Tue Jan 14 07:50:15 2014
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
  void operator()( tk::real* const particles, uint64_t npar, int nprop,
                   int offset, int ncomp ) {}
};

//! Zero initialization policy: zero particle properties
struct InitZero {
  void operator()( tk::real* const particles, uint64_t npar, int nprop,
                   int offset, int ncomp )
  {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (uint64_t p=0; p<npar; ++p) {
      memset( particles + p*nprop + offset, 0, ncomp*sizeof(tk::real) );
    }
  }
};

} // quinoa::

#endif // InitPolicy_h
