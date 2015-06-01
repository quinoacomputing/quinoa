//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 09:29:08 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Unit test driver
  \details   Unit test driver.
*/
//******************************************************************************
#ifndef UnitTestDriver_h
#define UnitTestDriver_h

#include "UnitTest/CmdLine/CmdLine.h"

namespace unittest {

class UnitTestPrint;

//! Unit test suite driver used polymorphically with tk::Driver
class UnitTestDriver {

  public:
    //! Constructor
    explicit UnitTestDriver( const UnitTestPrint& print,
                             const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}

  private:
    const UnitTestPrint& m_print;
};

} // unittest::

#endif // UnitTestDriver_h
