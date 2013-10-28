//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sun 27 Oct 2013 05:26:19 PM MDT
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
  public tk::Control< // tag    type
                      io,       ios > {

  public:
    //! Constructor: set all defaults
    CmdLine() {
      using namespace ctr;
      // Default I/O parameters
      set<io,control>("");
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
} // rngtest::

#endif // RNGTestCmdLine_h
