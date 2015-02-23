//******************************************************************************
/*!
  \file      src/Main/InciterDriver.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 09:04:59 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************
#ifndef InciterDriver_h
#define InciterDriver_h

#include <InciterPrint.h>
#include <Inciter/CmdLine/CmdLine.h>
#include <inciter.decl.h>

extern CProxy_Main mainProxy;

namespace inciter {

//! Inciter driver used polymorphically with tk::Driver
class InciterDriver {

  public:
    //! Constructor
    explicit InciterDriver( const InciterPrint& print,
                            const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() { mainProxy.finalize(); }

  private:
    const InciterPrint& m_print;
};

} // inciter::

#endif // InciterDriver_h
