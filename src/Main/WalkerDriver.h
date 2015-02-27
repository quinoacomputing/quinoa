//******************************************************************************
/*!
  \file      src/Main/WalkerDriver.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 09:46:18 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     WalkerDriver that drives Walker
  \details   WalkerDriver that drives Walker
*/
//******************************************************************************
#ifndef WalkerDriver_h
#define WalkerDriver_h

#include <WalkerPrint.h>
#include <Walker/CmdLine/CmdLine.h>

//! Everything that contributes to the walker executable
namespace walker {

//! Walker driver used polymorphically with Driver
class WalkerDriver {

  public:
    //! Constructor
    explicit WalkerDriver( const WalkerPrint& print,
                           const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}

  private:
    const WalkerPrint& m_print;
};

} // walker::

#endif // WalkerDriver_h
