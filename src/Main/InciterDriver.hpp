// *****************************************************************************
/*!
  \file      src/Main/InciterDriver.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter driver
  \details   Inciter driver.
*/
// *****************************************************************************
#ifndef InciterDriver_h
#define InciterDriver_h

#include "Inciter/CmdLine/CmdLine.hpp"

namespace inciter {

class InciterPrint;

//! Inciter driver used polymorphically with tk::Driver
class InciterDriver {

  public:
    //! Constructor
    explicit InciterDriver( const ctr::CmdLine& cmdline, int nrestart );

    //! Execute driver
    void execute() const;
};

} // inciter::

#endif // InciterDriver_h
