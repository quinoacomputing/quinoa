//******************************************************************************
/*!
  \file      src/Main/UnitTestDriver.h
  \author    J. Bakosi
  \date      Mon 11 Aug 2014 04:49:22 PM MDT
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
};

} // unittest::

#endif // UnitTestDriver_h
