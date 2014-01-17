//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Thu 16 Jan 2014 09:40:56 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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
      // Default I/O parameters
      set< tag::io, tag::control >( "" );
      set< tag::io, tag::input >( "" );
      set< tag::io, tag::output >( "out" );
      set< tag::io, tag::pdf >( "pdf" );
      set< tag::io, tag::glob >( "glob" );
      set< tag::io, tag::stat >( "stat" );
    }

    //! Destructor
    ~CmdLine() noexcept override = default;

    //! Instruct compiler to generate copy assigment
    CmdLine& operator=(const CmdLine&) = default;

  private:
    //! Don't permit copy constructor
    CmdLine(const CmdLine&) = delete;
    //! Don't permit move constructor
    CmdLine(CmdLine&&) = delete;
    //! Don't permit move assigment
    CmdLine& operator=(CmdLine&&) = delete;
};

//! CmdLine defaults
static const CmdLine CmdLineDefaults;

} // ctr::
} // quinoa::

#endif // QuinoaCmdLine_h
