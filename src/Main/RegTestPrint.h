//******************************************************************************
/*!
  \file      src/Main/RegTestPrint.h
  \author    J. Bakosi
  \date      Fri 20 Mar 2015 11:52:21 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     RegTest's pretty printer
  \details   RegTest's pretty printer
*/
//******************************************************************************
#ifndef RegTestPrint_h
#define RegTestPrint_h

#include <sstream>

#include <Types.h>
#include <Print.h>
#include <Exception.h>

namespace regtest {

//! RegTestPrint : tk::Print
class RegTestPrint : public tk::Print {

  public:
    //! Constructor
    //! \param[inout] str Verbose stream
    //! \param[inout] qstr Quiet stream
    //! \see tk::Print::Print
    //! \author J. Bakosi
    explicit RegTestPrint( std::ostream& str = std::clog,
                           std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}

};

} // regtest::

#endif // RegTestPrint_h
