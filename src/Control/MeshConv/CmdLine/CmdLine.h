//******************************************************************************
/*!
  \file      src/Control/MeshConv/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Mon 14 Jul 2014 09:30:38 PM MDT
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
} // meshconv::

#endif // MeshConvCmdLine_h
