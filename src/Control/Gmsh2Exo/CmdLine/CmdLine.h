//******************************************************************************
/*!
  \file      src/Control/Gmsh2Exo/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Wed Mar 19 10:14:30 2014
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Gmsh2Exo's command line
  \details   Gmsh2Exo's command line
*/
//******************************************************************************
#ifndef Gmsh2ExoCmdLine_h
#define Gmsh2ExoCmdLine_h

#include <string>

#include <Control.h>
#include <Gmsh2Exo/Types.h>

namespace gmsh2exo {
namespace ctr {

//! CmdLine : Control< specialized to Gmsh2Exo >, see Types.h,
class CmdLine :
  public tk::Control< // tag         type
                      tag::io,       ios > {

  public:
    //! Constructor: set all defaults
    CmdLine() {
      // Default I/O parameters
      set< tag::io, tag::control >("");
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
} // gmsh2exo::

#endif // Gmsh2ExoCmdLine_h
