//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 07:38:41 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
class CmdLine :
  public tk::Control< // tag         type
                      tag::io,       ios > {

  public:
    //! Constructor: set all defaults
    CmdLine() {
      // Default I/O parameters
      set< tag::io, tag::control >("");
    }
};

//! CmdLine defaults
static const CmdLine CmdLineDefaults;

} // ctr::
} // rngtest::

#endif // RNGTestCmdLine_h
