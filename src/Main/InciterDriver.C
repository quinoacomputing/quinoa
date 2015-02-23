//******************************************************************************
/*!
  \file      src/Main/InciterDriver.C
  \author    J. Bakosi
  \date      Mon 23 Feb 2015 08:32:21 AM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Inciter driver
  \details   Inciter driver.
*/
//******************************************************************************

#include <InciterDriver.h>
#include <Inciter/InputDeck/Parser.h>
//#include <integrator.decl.h>

namespace inciter {

extern ctr::InputDeck g_inputdeck;
//extern CProxy_Distributor g_DistributorProxy;

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
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );

  m_print.endpart();
  m_print.part( "Factory" );

  // Instantiate Distributor chare which drives the time-integration of
  // differential equations via several integrator chares. We only support a
  // single type of Distributor class at this point, so no factory
  // instantiation, simply fire up a Charm++ chare Distributor, which fires up
  // integrators. Store proxy handle in global-scope to make it available to
  // individual integrators so they can call back to Distributor. Since this
  // is called inside the main chare constructor, the Charm++ runtime system
  // distributes the handle along with all other global-scope data.
  //g_DistributorProxy = CProxy_Distributor::ckNew( cmdline );
}
