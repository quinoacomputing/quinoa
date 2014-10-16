//******************************************************************************
/*!
  \file      src/Main/QuinoaDriver.C
  \author    J. Bakosi
  \date      Tue 26 Aug 2014 12:30:20 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     QuinoaDriver that drives Quinoa
  \details   QuinoaDriver that drives Quinoa
*/
//******************************************************************************

#include <QuinoaDriver.h>
#include <Quinoa/InputDeck/Parser.h>
#include <integrator.decl.h>

namespace quinoa {

extern ctr::InputDeck g_inputdeck;
extern CProxy_Distributor g_DistributorProxy;

} // quinoa::

using quinoa::QuinoaDriver;

QuinoaDriver::QuinoaDriver( const QuinoaPrint& print,
                            const ctr::CmdLine& cmdline ) :
  m_print( print )
//******************************************************************************
//  Constructor
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
  g_DistributorProxy = CProxy_Distributor::ckNew( cmdline );
}
