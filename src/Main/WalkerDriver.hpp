// *****************************************************************************
/*!
  \file      src/Main/WalkerDriver.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     WalkerDriver that drives Walker
  \details   WalkerDriver that drives Walker
*/
// *****************************************************************************
#ifndef WalkerDriver_h
#define WalkerDriver_h

#include "Walker/CmdLine/CmdLine.hpp"

//! Everything that contributes to the walker executable
namespace walker {

//! Walker driver used polymorphically with Driver
class WalkerDriver {

  public:
    //! Constructor
    explicit WalkerDriver( const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}
};

} // walker::

#endif // WalkerDriver_h
