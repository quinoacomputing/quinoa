// *****************************************************************************
/*!
  \file      src/Main/WalkerDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     WalkerDriver that drives Walker
  \details   WalkerDriver that drives Walker
*/
// *****************************************************************************

#include "Tags.hpp"
#include "WalkerPrint.hpp"
#include "WalkerDriver.hpp"
#include "Walker/InputDeck/Parser.hpp"
#include "Walker/CmdLine/CmdLine.hpp"
#include "Walker/InputDeck/InputDeck.hpp"
#include "TaggedTupleDeepPrint.hpp"
#include "Writer.hpp"

#include "NoWarning/distributor.decl.h"

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern CProxy_Distributor g_DistributorProxy;

} // walker::

using walker::WalkerDriver;

WalkerDriver::WalkerDriver( const ctr::CmdLine& cmdline, int nrestart )
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
    g_inputdeck_defaults.get< tag::cmd, tag::io, tag::screen >();
  WalkerPrint print( cmdline.logname( def, nrestart ),
                     cmdline.get< tag::verbose >() ? std::cout : std::clog,
                     std::ios_base::app );

  // Parse input deck into g_inputdeck
  print.item( "Control file", cmdline.get< tag::io, tag::control >() );
  InputDeckParser inputdeckParser( print, cmdline, g_inputdeck );
  print.item( "Parsed control file", "success" );
  print.endpart();

  // Output command line object to file
  auto logfilename = tk::walker_executable() + "_input.log";
  tk::Writer log( logfilename );
  tk::print( log.stream(), "inputdeck", g_inputdeck );

  // Instantiate Distributor chare on PE 0 which drives the time-integration of
  // differential equations via several integrator chares. We only support a
  // single type of Distributor class at this point, so no factory
  // instantiation, simply fire up a Charm++ chare Distributor, which fires up
  // integrators. Store proxy handle in global-scope to make it available to
  // individual integrators so they can call back to Distributor. Since this
  // is called inside the main chare constructor, the Charm++ runtime system
  // distributes the handle along with all other global-scope data.
  g_DistributorProxy = CProxy_Distributor::ckNew( 0 );
}
