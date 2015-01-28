//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 11:43:15 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Unit test driver
  \details   Unit test driver.
*/
//******************************************************************************
#ifndef UnitTestDriver_h
#define UnitTestDriver_h

#include <UnitTestPrint.h>
#include <UnitTest/CmdLine/CmdLine.h>

namespace unittest {

//! Unit test suite driver used polymorphically with tk::Driver
class UnitTestDriver {

  public:
    //! Constructor
    explicit UnitTestDriver( const UnitTestPrint& print,
                             const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() {}

  private:
    const UnitTestPrint& m_print;
};

} // unittest::

#endif // UnitTestDriver_h
