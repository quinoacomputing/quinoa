//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.h
  \author    J. Bakosi
  \date      Thu 19 Nov 2015 09:48:32 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Functions used to setup inciter
  \details   Functions used to setup inciter.
*/
//******************************************************************************
#ifndef InciterSetup_h
#define InciterSetup_h

#include <vector>
#include <map>
#include <cstddef>
#include <iosfwd>
#include <utility>

#include "Timer.h"
#include "InciterPrint.h"
#include "Inciter/InputDeck/InputDeck.h"
#include "Inciter/CmdLine/CmdLine.h"

namespace tk { class Print; }

namespace inciter {

// Parse command line
void
parseCmdLine( int argc, char** argv, ctr::CmdLine& cmdline );

//! Parse command line and input deck, instantiate pretty printer, echo info
void
init( const ctr::CmdLine& cmdline,
      const tk::Print& print,
      ctr::InputDeck& inputdeck,
      int argc,
      char** argv );

//! Prepare computational mesh
std::map< int, std::vector< std::size_t > >
prepareMesh(
  const ctr::InputDeck& inputdeck,
  const InciterPrint& print,
  std::vector< std::pair< std::string, tk::Timer::Watch > >& timestamp,
  uint64_t& nchare,
  std::size_t& npoin );

} // inciter::

#endif // InciterSetup_h
