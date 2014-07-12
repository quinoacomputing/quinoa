//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed 11 Jun 2014 01:55:33 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     MeshConv's command line
  \details   MeshConv's command line
*/
//******************************************************************************
#ifndef MeshConvCmdLine_h
#define MeshConvCmdLine_h

#include <string>

#include <Control.h>
#include <MeshConv/Types.h>

namespace meshconv {
namespace ctr {

//! CmdLine : Control< specialized to MeshConv >, see Types.h,
class CmdLine :
  public tk::Control< // tag         type
                      tag::io,       ios > {
};

} // ctr::
} // meshconv::

#endif // MeshConvCmdLine_h
