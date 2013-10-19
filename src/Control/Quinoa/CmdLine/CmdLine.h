//******************************************************************************
/*!
  \file      src/Control/Quinoa/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Sat 19 Oct 2013 08:01:29 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's command line
  \details   Quinoa's command line
*/
//******************************************************************************
#ifndef QuinoaCmdLine_h
#define QuinoaCmdLine_h

#include <string>

#include <Control.h>

namespace quinoa {
namespace ctr {

//! CmdLine : Control< specialized to Quinoa >, see Types.h,
class CmdLine :
  public tk::Control< // tag    type
                      io,       ios > {

  public:
    //! Constructor: set all defaults
    CmdLine() {
      using namespace ctr;
      // Default I/O parameters
      set<io,control>("");
      set<io,input>("");
      set<io,output>("out");
      set<io,pdf>("pdf");
      set<io,glob>("glob");
      set<io,stat>("stat");
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
