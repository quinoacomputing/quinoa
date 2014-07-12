//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 10:12:33 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Driver random number test suite driver
  \details   Driver random number test suite driver
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
