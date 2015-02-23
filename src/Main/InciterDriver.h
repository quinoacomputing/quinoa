//******************************************************************************
/*!
  \file      src/Main/InciterDriver.h
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 03:18:29 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************
#ifndef InciterDriver_h
#define InciterDriver_h

#include <InciterPrint.h>
#include <Inciter/CmdLine/CmdLine.h>
#include <UnsMesh.h>

namespace inciter {

//! Inciter driver used polymorphically with tk::Driver
class InciterDriver {

  public:
    //! Constructor
    explicit InciterDriver( const InciterPrint& print,
                            const ctr::CmdLine& cmdline );

    //! Execute driver
    void execute() const;

  private:
    //! Echo mesh statistics
    void meshStats( const tk::UnsMesh& mesh ) const;

    const InciterPrint& m_print;        //!< Pretty printer
    std::string m_input;                //!< Input file name
};

} // inciter::

#endif // InciterDriver_h
