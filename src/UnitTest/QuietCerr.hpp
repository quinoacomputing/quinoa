// *****************************************************************************
/*!
  \file      src/UnitTest/QuietCerr.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ nodegroup to quiet std::cerr in a thread-safe fashion
  \details   Charm++ nodegroup to quiet std::cerr in a thread-safe fashion.
*/
// *****************************************************************************
#ifndef QuietCerr_h
#define QuietCerr_h

#include "NoWarning/quietcerr.decl.h"

namespace tk {

//! Chare state Charm++ chare nodegroup class
//! \details Instantiations of QuietCerr comprise a processor aware Charm++
//!   chare node group. When instantiated, a new object is created on each
//!   compute node and not more (as opposed to individual chares or chare array
//!   object elements). See also the Charm++ interface file quietcerr.ci.
class QuietCerr : public CBase_QuietCerr {

  public:
    //! Quiet std::cerr by redirecting its stream state to a stringstream
    static void quiet();

    //! Destructor: restore std::cerr's stream state
    ~QuietCerr() override;
};

} // tk::

#endif // QuietCerr_h
