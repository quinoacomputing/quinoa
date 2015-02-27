//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 11:42:12 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Random number generator test suite driver
  \details   Random number generator test suite driver.
*/
//******************************************************************************
#ifndef RNGTestDriver_h
#define RNGTestDriver_h

#include <Battery.h>
#include <RNGTestPrint.h>

//! Everything that contributes to the rngtest executable
namespace rngtest {

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
