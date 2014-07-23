//******************************************************************************
/*!
  \file      src/Main/UnitTestPrint.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 09:58:19 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTest's printer
  \details   UnitTest's printer
*/
//******************************************************************************
#ifndef UnitTestPrint_h
#define UnitTestPrint_h

#include <Types.h>
#include <Print.h>
#include <flip_map.h>

namespace unittest {

//! UnitTestPrint : Print
class UnitTestPrint : public tk::Print {

  public:
    //! Constructor
    explicit UnitTestPrint( std::ostream& str = tk::null,
                            std::ostream& qstr = std::cout ) :
      Print( str, qstr ) {}
};

} // unittest::

#endif // UnitTestPrint_h
