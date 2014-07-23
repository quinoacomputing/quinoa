//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:20:13 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     UnitTest's command line
  \details   UnitTest's command line
*/
//******************************************************************************
#ifndef UnitTestCmdLine_h
#define UnitTestCmdLine_h

#include <string>

#include <Control.h>
#include <UnitTest/Types.h>

namespace unittest {
namespace ctr {

//! CmdLine : Control< specialized to UnitTest >, see Types.h,
class CmdLine :
  public tk::Control< // tag            type
                      tag::io,          ios,
                      tk::tag::verbose, bool > {
  public:
    //! Constructor: set defaults. Anything not set here initialized by the
    //! compiler using the default constructor for the corresponding type.
    CmdLine() {
      set< tk::tag::verbose >( false ); // Use quiet output by default
    }
};

} // ctr::
} // unittest::

#endif // UnitTestCmdLine_h
