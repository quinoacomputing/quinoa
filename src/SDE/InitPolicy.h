//******************************************************************************
/*!
  \file      src/SDE/InitPolicy.h
  \author    J. Bakosi
  \date      Thu 06 Feb 2014 05:18:52 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Initialization policies
  \details   Initialization policies
*/
//******************************************************************************
#ifndef InitPolicy_h
#define InitPolicy_h

#include <cstring>

#include <boost/mpl/vector.hpp>

#include <Types.h>
#include <Quinoa/Options/InitPolicy.h>

namespace quinoa {

//! Raw initialization policy: leave memory uninitialized
struct InitRaw {

  InitRaw() = default;

  std::string policy() const noexcept {
    return ctr::InitPolicy().name( ctr::InitPolicyType::RAW );
  }

  ctr::InitPolicyType type() const noexcept {
    return ctr::InitPolicyType::RAW;
  }

  InitRaw( const ParProps& particles, uint64_t npar, uint32_t nprop, int offset,
           int ncomp, int nthreads ) {}
};

//! Zero initialization policy: zero particle properties
struct InitZero {

  InitZero() = default;

  std::string policy() const noexcept {
    return ctr::InitPolicy().name( ctr::InitPolicyType::ZERO );
  }

  ctr::InitPolicyType type() const noexcept {
    return ctr::InitPolicyType::ZERO;
  }

  InitZero( const ParProps& particles, uint64_t npar, uint32_t nprop,
            int offset, int ncomp, int nthreads )
  {
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

//! List of all initialization policies
using InitPolicies = boost::mpl::vector< InitRaw, InitZero >;

} // quinoa::

#endif // InitPolicy_h
