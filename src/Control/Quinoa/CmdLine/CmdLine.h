//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed 06 Aug 2014 08:51:20 PM MDT
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
class CmdLine : public tk::Control< // tag            type
                                    tag::io,          ios,
                                    tk::tag::verbose, bool > {

  public:
    //! Constructor: set all defaults. Anything not set here initialized by the
    //! compiler using the default constructor for the corresponding type.
    CmdLine() {
      set< tag::io, tag::output >( "out" );
      set< tag::io, tag::pdf >( "pdf" );
      set< tag::io, tag::glob >( "glob" );
      set< tag::io, tag::stat >( "stat" );
      set< tk::tag::verbose >( false ); // Use quiet output by default
    }
};

} // ctr::
} // quinoa::

#endif // QuinoaCmdLine_h
