// *****************************************************************************
/*!
  \file      src/Main/UnitTestDriver.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Unit test driver
  \details   Unit test driver.
*/
// *****************************************************************************
#ifndef UnitTestDriver_h
#define UnitTestDriver_h

#include "UnitTest/CmdLine/CmdLine.hpp"

namespace unittest {

//! Unit test suite driver used polymorphically with tk::Driver
class UnitTestDriver {

  public:
    //! Constructor
    explicit UnitTestDriver( const ctr::CmdLine& cmdline, int );

    //! Execute driver
    void execute() const {}
};

} // unittest::

#endif // UnitTestDriver_h
