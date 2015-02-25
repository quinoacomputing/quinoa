//******************************************************************************
/*!
  \file      src/Main/InciterDriver.h
  \author    J. Bakosi
  \date      Tue 24 Feb 2015 12:16:36 PM MST
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
    //! Print information at startup
    void info( const tk::UnsMesh& mesh,
               uint64_t chunksize,
               uint64_t remainder,
               uint64_t nchare ) const;

    //! Partition mesh using Zoltan
    void partition( const tk::UnsMesh& mesh ) const;

    const InciterPrint& m_print;        //!< Pretty printer
    std::string m_input;                //!< Input file name
};

} // inciter::

#endif // InciterDriver_h
