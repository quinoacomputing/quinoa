//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed 11 Jun 2014 01:54:50 PM MDT
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
class CmdLine :
  public tk::Control< // tag         type
                      tag::io,       ios > {
};

} // ctr::
} // rngtest::

#endif // RNGTestCmdLine_h
