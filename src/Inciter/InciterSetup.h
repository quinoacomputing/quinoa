//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.h
  \author    J. Bakosi
  \date      Sat 21 Nov 2015 02:58:22 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Functions used to setup inciter
  \details   Functions used to setup inciter.
*/
//******************************************************************************
#ifndef InciterSetup_h
#define InciterSetup_h

#include <vector>
#include <unordered_map>
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
std::unordered_map< int, std::vector< std::size_t > >
prepareMesh(
  const ctr::InputDeck& inputdeck,
  const InciterPrint& print,
  std::vector< std::pair< std::string, tk::Timer::Watch > >& timestamp,
  uint64_t& nchare,
  std::size_t& npoin );

} // inciter::

#endif // InciterSetup_h
