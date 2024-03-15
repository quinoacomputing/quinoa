// *****************************************************************************
/*!
  \file      src/Main/InciterDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter driver
  \details   Inciter driver.
*/
// *****************************************************************************

#include "InciterPrint.hpp"
#include "InciterDriver.hpp"
#include "Inciter/InputDeck/Parser.hpp"
#include "Inciter/InputDeck/LuaParser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"

#include "NoWarning/transporter.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;
extern ctr::InputDeck g_oldinputdeck_defaults;
extern ctr::New2InputDeck g_inputdeck_defaults;
extern ctr::New2InputDeck g_newinputdeck;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const ctr::CmdLine& cmdline, int nrestart )
// *****************************************************************************
//  Constructor
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \param[in] nrestart Number of times restarted
// *****************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

  // Create pretty printer
  const auto& def =
    g_inputdeck_defaults.get< newtag::cmd, tag::io, tag::screen >();
  InciterPrint print( cmdline.logname( def, nrestart ),
                      cmdline.get< tag::verbose >() ? std::cout : std::clog,
                      std::ios_base::app );

  print.item( "Non-blocking migration, -" + *kw::nonblocking::alias(),
               cmdline.get< tag::nonblocking >() ? "on" : "off" );
  print.item( "Benchmark mode, -" + *kw::benchmark::alias(),
               cmdline.get< tag::benchmark >() ? "on" : "off" );

  auto lbfreq = cmdline.get< tag::lbfreq >();
  if ( lbfreq < kw::lbfreq::info::expect::lower ||
       lbfreq > kw::lbfreq::info::expect::upper ) {
    Throw( "Load-balancing frequency should be greater than 0." );
  }
  print.item( "Load-balancing frequency, -" + *kw::lbfreq::alias(),
               std::to_string(cmdline.get< tag::lbfreq >()) );

  auto rsfreq = cmdline.get< tag::rsfreq >();
  if ( rsfreq < kw::rsfreq::info::expect::lower ||
       rsfreq > kw::rsfreq::info::expect::upper ) {
    Throw( "Checkpoint/restart frequency should be greater than 0." );
  }
  print.item( "Checkpoint/restart frequency, -" + *kw::rsfreq::alias(),
               std::to_string(cmdline.get< tag::rsfreq >()) );

  // Parse input deck into g_inputdeck
  print.item( "Control file", cmdline.get< tag::io, tag::control >() );
  g_newinputdeck = g_inputdeck_defaults;   // overwrite with defaults if restarted

  LuaParser luaparser( print, cmdline, g_newinputdeck );
  print.item( "Parsed lua file", "success" );

  //InputDeckParser inputdeckParser( print, cmdline, g_inputdeck );
  //print.item( "Parsed control file", "success" );
  print.endpart();

  // Output command line object to file
  auto logfilename = tk::inciter_executable() + "_input.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), "inputdeck", g_inputdeck );

  CkExit();
}

void
InciterDriver::execute() const
// *****************************************************************************
//  Run inciter
// *****************************************************************************
{
  // Instantiate Transporter chare on PE 0 which drives time-integration
  CProxy_Transporter::ckNew( 0 );
}
