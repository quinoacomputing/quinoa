//******************************************************************************
/*!
  \file      src/Inciter/InciterSetup.h
  \author    J. Bakosi
  \date      Sat 30 May 2015 10:45:37 AM MDT
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
ctr::CmdLine
parseCmdLine( int argc, char** argv );

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
  std::size_t& npoin,
  std::vector< std::vector< std::size_t > >& point,
  std::vector< std::vector< std::size_t > >& element,
  std::vector< std::size_t >& meshfilemap,
  std::vector< std::size_t >& tetinpoel,
  std::pair< std::vector< std::size_t >, std::vector< std::size_t > >& esup,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& pcomm,
  std::vector< std::map< std::size_t, std::vector< std::size_t > > >& ecomm );

} // inciter::

#endif // InciterSetup_h
