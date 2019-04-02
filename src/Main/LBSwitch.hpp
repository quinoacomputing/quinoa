// *****************************************************************************
/*!
  \file      src/Main/LBSwitch.hpppp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Charm++ chare group for switching on/off load balancing
  \details   Charm++ chare group for switching on/off load balancing.
*/
// *****************************************************************************
#ifndef LBSwitch_h
#define LBSwitch_h

#include "NoWarning/lbswitch.decl.h"

namespace tk {

//! Load balancer switch Charm++ chare group class
//! \details Instantiations of LBSwitch comprise a processor aware
//!   Charm++ chare group. When instantiated, a new object is created on each
//!   PE and not more (as opposed to individual chares or chare array object
//!   elements). See also the Charm++ interface file lbswitch.ci.
class LBSwitch : public CBase_LBSwitch {

  public:
    //! Constructor: turn on automatic load balancing
    explicit LBSwitch( bool verbose );

    //! Turn off automatic load balancing
    static void off();
};

} // tk::

#endif // LBSwitch_h
