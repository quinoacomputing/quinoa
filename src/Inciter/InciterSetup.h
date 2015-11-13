//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.h
  \author    J. Bakosi
  \date      Mon 09 Nov 2015 12:31:38 PM MST
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

namespace tk { class Print; }

namespace inciter {

namespace ctr { class InputDeck; class CmdLine; }

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
void
prepareMesh(
  const ctr::CmdLine& cmdline,
  const tk::Print& print,
  const ctr::InputDeck& inputdeck,
  std::vector< std::pair< std::string, tk::Timer::Watch > >& timestamp,
  uint64_t& nchare,
  std::size_t& npoin,
  std::map< int, std::vector< std::size_t > >& element );

} // inciter::

#endif // InciterSetup_h
