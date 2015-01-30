//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Wed 28 Jan 2015 11:40:15 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Quinoa driver
  \details   Quinoa driver.
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <QuinoaPrint.h>
#include <Quinoa/CmdLine/CmdLine.h>

namespace quinoa {

//! Quinoa driver used polymorphically with tk::Driver
class QuinoaDriver {

  public:
    //! Constructor
    explicit QuinoaDriver( const QuinoaPrint& print,
                           const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() {}

  private:
    const QuinoaPrint& m_print;
};

} // quinoa::

#endif // QuinoaDriver_h
