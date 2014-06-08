//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sun 08 Jun 2014 03:42:24 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

  public:
    //! Constructor: set all defaults
    CmdLine() {
      // Default I/O parameters
      set< tag::io, tag::input >( "" );
      set< tag::io, tag::output >( "" );
    }
};

//! CmdLine defaults
static const CmdLine CmdLineDefaults;

} // ctr::
} // meshconv::

#endif // MeshConvCmdLine_h
