//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed 11 Jun 2014 01:51:25 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's command line
  \details   Quinoa's command line
*/
//******************************************************************************
#ifndef QuinoaCmdLine_h
#define QuinoaCmdLine_h

#include <string>

#include <Control.h>
#include <Quinoa/Types.h>

namespace quinoa {
namespace ctr {

//! CmdLine : Control< specialized to Quinoa >, see Types.h,
class CmdLine : public tk::Control< // tag    type
                                    tag::io,  ios > {

  public:
    //! Constructor: set all defaults
    CmdLine() {
      set< tag::io, tag::control >( "" );
      set< tag::io, tag::input >( "" );
      set< tag::io, tag::output >( "out" );
      set< tag::io, tag::pdf >( "pdf" );
      set< tag::io, tag::glob >( "glob" );
      set< tag::io, tag::stat >( "stat" );
    }
};

} // ctr::
} // quinoa::

#endif // QuinoaCmdLine_h
