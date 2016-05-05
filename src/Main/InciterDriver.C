//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Tue 03 May 2016 09:41:26 AM MDT
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************

#include <unordered_map>

#include "InciterPrint.h"
#include "InciterDriver.h"
#include "Inciter/InputDeck/Parser.h"
#include "Inciter/CmdLine/CmdLine.h"
#include "Inciter/InputDeck/InputDeck.h"

#include "NoWarning/conductor.decl.h"

namespace inciter {

extern ctr::InputDeck g_inputdeck;

} // inciter::

using inciter::InciterDriver;

InciterDriver::InciterDriver( const InciterPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
//******************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here (if any)

  // Parse input deck into g_inputdeck
  m_print.item( "Control file", cmdline.get< tag::io, tag::control >() );  
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );
  m_print.item( "Parsed control file", "success" );  
  m_print.endpart();
}

void
InciterDriver::execute() const
//******************************************************************************
//  Run inciter
//! \author J. Bakosi
//******************************************************************************
{
  // Instantiate Conductor chare which drives the time-integration of a PDE via
  // several Performer chares.
  CProxy_Conductor::ckNew();
}
