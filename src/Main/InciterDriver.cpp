// *****************************************************************************
/*!
  \file      src/Main/InciterDriver.cpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Inciter driver
  \details   Inciter driver.
*/
// *****************************************************************************

#include <unordered_map>

#include "InciterPrint.hpp"
#include "InciterDriver.hpp"
#include "Inciter/InputDeck/Parser.hpp"
#include "Inciter/CmdLine/CmdLine.hpp"
#include "Inciter/InputDeck/InputDeck.hpp"

#include "NoWarning/transporter.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const InciterPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print )
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
// *****************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

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

  // Parse input deck into g_inputdeck
  m_print.item( "Control file", cmdline.get< tag::io, tag::control >() );  
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );
  m_print.item( "Parsed control file", "success" );  
  m_print.endpart();
}

void
InciterDriver::execute() const
// *****************************************************************************
//  Run inciter
// *****************************************************************************
{
  // Instantiate Transporter chare on PE 0 which drives the time-integration of
  // a PDE via several worker chares.
  CProxy_Transporter::ckNew( 0 );
}
