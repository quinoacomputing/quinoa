//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Fri 08 Jan 2016 06:11:29 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
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

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "conductor.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

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
