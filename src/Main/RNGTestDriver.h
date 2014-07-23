//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Tue 22 Jul 2014 10:21:39 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator test suite driver
  \details   Random number generator test suite driver
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

//! Random number generator test suite driver used polymorphically with Driver
class RNGTestDriver {

  public:
    //! Constructor
    explicit RNGTestDriver( const RNGTestPrint& print,
                            const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute();

  private:
    const RNGTestPrint& m_print;
};

} // rngtest::

#endif // RNGTestDriver_h
