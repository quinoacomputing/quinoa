//******************************************************************************
/*!
  \file      src/Base/LoadDistributor.h
  \author    J. Bakosi
  \date      Tue 24 Feb 2015 10:14:22 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Load distributors
  \details   Load distributors compute chunksize based on the degree of
     virtualization.
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
                       uint64_t& chunksize,
                       uint64_t& remainder );

} // tk::

#endif // LoadDistributor_h
