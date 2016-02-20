//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 09:15:01 AM MDT
  \copyright 2012-2016, Jozsef Bakosi.
  \brief     Random number generator test suite driver
  \details   Random number generator test suite driver.
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <map>
#include <functional>

#include "RNGTest/Options/Battery.h"
#include "RNGTest/CmdLine/CmdLine.h"

//! Everything that contributes to the rngtest executable
namespace rngtest {

class Battery;
class RNGTestPrint;

//! Battery factory type
using BatteryFactory = std::map< ctr::BatteryType, std::function< Battery() > >;

//! \brief Random number generator test suite driver used polymorphically with
//!   tk::Driver
class RNGTestDriver {

  public:
    //! Constructor
    explicit RNGTestDriver( const RNGTestPrint& print,
                            const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const;

  private:
    const RNGTestPrint& m_print;
};

} // rngtest::

#endif // RNGTestDriver_h
