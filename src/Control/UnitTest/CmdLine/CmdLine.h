//******************************************************************************
/*!
  \file      src/Control/UnitTest/CmdLine/CmdLine.h
  \author    J. Bakosi
  \date      Mon 08 Dec 2014 02:32:43 PM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     UnitTest's command line
  \details   UnitTest's command line
*/
//******************************************************************************
#ifndef UnitTestCmdLine_h
#define UnitTestCmdLine_h

#include <string>

#include <Control.h>
#include <UnitTest/Types.h>

namespace unittest {
namespace ctr {

//! CmdLine : Control< specialized to UnitTest >, see Types.h,
class CmdLine : public tk::Control<
                  // tag        type
                  tag::verbose, bool,
                  tag::error,   std::vector< std::string > > {
  public:
    //! Constructor: set defaults. Anything not set here initialized by the
    //! compiler using the default constructor for the corresponding type.
    CmdLine() {
      set< tag::verbose >( false ); // Use quiet output by default
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::verbose, bool,
                   tag::error,   std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, CmdLine& c ) { c.pup(p); }
};

} // ctr::
} // unittest::

#endif // UnitTestCmdLine_h
