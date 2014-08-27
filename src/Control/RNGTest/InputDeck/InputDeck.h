//******************************************************************************
/*!
  \file      src/Control/RNGTest/InputDeck/InputDeck.h
  \author    J. Bakosi
  \date      Fri 22 Aug 2014 10:53:54 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Random number generator test suite input deck
  \details   Random number generator test suite input deck
*/
//******************************************************************************
#ifndef RNGTestInputDeck_h
#define RNGTestInputDeck_h

#include <Control.h>
#include <RNGTest/CmdLine/CmdLine.h>

namespace rngtest {
namespace ctr {

//! InputDeck : Control< specialized to RNGTest >, see Types.h
//! This is also where the command line parser stores
class InputDeck : public tk::Control<
                    // tag           type
                    tag::title,      std::string,
                    tag::selected,   selects,
                    tag::io,         ios,
                    tag::cmd,        CmdLine,
                    tag::param,      parameters,
                    tk::tag::error,  std::vector< std::string > > {

  public:
    //! Pack/Unpack
    void pup( PUP::er& p ) {
      tk::Control< tag::title,      std::string,
                   tag::selected,   selects,
                   tag::io,         ios,
                   tag::cmd,        CmdLine,
                   tag::param,      parameters,
                   tk::tag::error,  std::vector< std::string > >::pup(p);
    }
    friend void operator|( PUP::er& p, InputDeck& c ) { c.pup(p); }
};

} // ctr::
} // rngtest::

#endif // RNGTestInputDeck_h
