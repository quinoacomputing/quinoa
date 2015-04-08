//******************************************************************************
/*!
  \file      src/Base/LoadDistributor.h
  \author    J. Bakosi
  \date      Mon 06 Apr 2015 11:04:31 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Load distributors and partitioning data types
  \details   Load distributors and partitioning data types. Load distributors
     compute chunksize based on the degree of virtualization.
*/
//******************************************************************************
#ifndef LoadDistributor_h
#define LoadDistributor_h

#include <sstream>

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
