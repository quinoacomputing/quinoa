//******************************************************************************
/*!
  \file      src/Main/InciterDriver.h
  \author    J. Bakosi
  \date      Thu 05 Nov 2015 02:59:57 PM MST
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
    void execute(
      std::size_t npoin,
      uint64_t nchare,
      const std::vector< std::vector< std::vector< std::size_t > > >& gelement )
    const;

  private:
    const InciterPrint& m_print;        //!< Pretty printer
};

} // inciter::

#endif // InciterDriver_h
