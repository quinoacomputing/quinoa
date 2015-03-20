//******************************************************************************
/*!
  \file      src/Main/RegTestDriver.h
  \author    J. Bakosi
  \date      Fri 20 Mar 2015 11:31:21 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Regression test harness driver
  \details   Regression test harness driver.
*/
//******************************************************************************
#ifndef RegTestDriver_h
#define RegTestDriver_h

#include <RegTestPrint.h>
#include <RegTest/CmdLine/CmdLine.h>

namespace regtest {

//! Regression test suite driver used polymorphically with tk::Driver
class RegTestDriver {

  public:
    //! Constructor
    explicit RegTestDriver( const RegTestPrint& print,
                            const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}

  private:
    const RegTestPrint& m_print;
};

} // regtest::

#endif // RegTestDriver_h
