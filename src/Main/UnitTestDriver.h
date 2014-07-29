//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.h
  \author    J. Bakosi
  \date      Tue 29 Jul 2014 08:05:41 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Unit test driver
  \details   Unit test driver
*/
//******************************************************************************
#ifndef UnitTestDriver_h
#define UnitTestDriver_h

#include <UnitTestPrint.h>
#include <UnitTest/CmdLine/CmdLine.h>

//! Everything that contributes to the unittest executable
namespace unittest {

//! Unit test suite driver used polymorphically with Driver
class UnitTestDriver {

  public:
    //! Constructor
    explicit UnitTestDriver( const UnitTestPrint& print,
                             const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() {}

  private:
    const UnitTestPrint& m_print;
    const ctr::CmdLine& m_cmdline;
};

} // unittest::

#endif // UnitTestDriver_h
