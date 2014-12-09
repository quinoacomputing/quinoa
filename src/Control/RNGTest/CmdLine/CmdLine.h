//******************************************************************************
/*!
  \file      src/Control/RNGTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:25:45 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
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
class CmdLine : public tk::Control<
                  // tag        type
                  tag::io,      ios,
                  tag::verbose, bool,
                  tag::error,   std::vector< std::string > > {
  public:
    //! Constructor: set defaults. Anything not set here initialized by the
    //! compiler using the default constructor for the corresponding type.
    CmdLine() {
      set< tag::verbose >( false ); // Use quiet output by default
    }
};

} // ctr::
} // rngtest::

#endif // RNGTestCmdLine_h
