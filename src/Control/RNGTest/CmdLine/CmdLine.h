//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sat 12 Jul 2014 10:43:33 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNGTest's command line
  \details   RNGTest's command line
*/
//******************************************************************************
#ifndef RNGTestCmdLine_h
#define RNGTestCmdLine_h

#include <string>

#include <Control.h>
#include <RNGTest/Types.h>

namespace rngtest {
namespace ctr {

//! CmdLine : Control< specialized to RNGTest >, see Types.h,
class CmdLine : public tk::Control< // tag            type
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
} // rngtest::

#endif // RNGTestCmdLine_h
