//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.h
  \author    J. Bakosi
  \date      Fri 15 Aug 2014 10:12:21 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************
#ifndef QuinoaDriver_h
#define QuinoaDriver_h

#include <QuinoaPrint.h>
#include <Quinoa/CmdLine/CmdLine.h>

//! Everything that contributes to the quinoa executable
namespace quinoa {

//! Quinoa driver used polymorphically with Driver
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
