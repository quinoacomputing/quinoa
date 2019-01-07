// *****************************************************************************
/*!
  \file      src/Main/UnitTestDriver.h
  \copyright 2012-2015, J. Bakosi, 2016-2018, Los Alamos National Security, LLC.
  \brief     Unit test driver
  \details   Unit test driver.
*/
// *****************************************************************************
#ifndef UnitTestDriver_h
#define UnitTestDriver_h

#include "UnitTest/CmdLine/CmdLine.h"

namespace unittest {

class UnitTestPrint;

//! Unit test suite driver used polymorphically with tk::Driver
class UnitTestDriver {

  public:
    //! Constructor
    explicit UnitTestDriver( const UnitTestPrint&,
                             const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}
};

} // unittest::

#endif // UnitTestDriver_h
