//******************************************************************************
/*!
  \file      src/Main/InciterDriver.h
  \author    J. Bakosi
  \date      Sun 31 May 2015 06:08:35 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************
#ifndef InciterDriver_h
#define InciterDriver_h

namespace inciter {

class InciterPrint;

//! Inciter driver used polymorphically with tk::Driver
class InciterDriver {

  public:
    //! Constructor
    explicit InciterDriver( const InciterPrint& print );

    //! Execute driver
    void execute() const;

  private:
    const InciterPrint& m_print;        //!< Pretty printer
};

} // inciter::

#endif // InciterDriver_h
