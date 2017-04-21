// *****************************************************************************
/*!
  \file      src/LoadBalance/LoadDistributor.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Load distributors and partitioning data types
  \details   Load distributors and partitioning data types. Load distributors
     compute chunksize based on the degree of virtualization.
*/
// *****************************************************************************
#ifndef LoadDistributor_h
#define LoadDistributor_h

#include <cstdint>

#include "Types.h"

namespace tk {

//! Compute linear load distribution for given total work and virtualization
uint64_t
linearLoadDistributor( tk::real virtualization,
                       uint64_t load,
                       int npe,
                       uint64_t& chunksize,
                       uint64_t& remainder );

} // tk::

#endif // LoadDistributor_h
