// *****************************************************************************
/*!
  \file      src/Base/LoadDistributor.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Load distributors and partitioning data types
  \details   Load distributors and partitioning data types. Load distributors
     compute chunksize based on the degree of virtualization.
*/
// *****************************************************************************
#ifndef LoadDistributor_h
#define LoadDistributor_h

#include <cstdint>

#include "Types.hpp"

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
