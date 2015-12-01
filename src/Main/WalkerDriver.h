//******************************************************************************
/*!
  \file      src/Main/WalkerDriver.h
  \author    J. Bakosi
  \date      Mon 30 Nov 2015 12:34:32 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     WalkerDriver that drives Walker
  \details   WalkerDriver that drives Walker
*/
//******************************************************************************
#ifndef WalkerDriver_h
#define WalkerDriver_h

#include "Walker/CmdLine/CmdLine.h"

//! Everything that contributes to the walker executable
namespace walker {

class WalkerPrint;

//! Walker driver used polymorphically with Driver
class WalkerDriver {

  public:
    //! Constructor
    explicit WalkerDriver( const WalkerPrint& print,
                           const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const {}

  private:
    const WalkerPrint& m_print;        //!< Pretty printer
};

} // walker::

#endif // WalkerDriver_h
