//******************************************************************************
/*!
  \file      src/Main/WalkerDriver.C
  \author    J. Bakosi
  \date      Mon 01 Jun 2015 09:43:36 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     WalkerDriver that drives Walker
  \details   WalkerDriver that drives Walker
*/
//******************************************************************************

#include <string>

#include <boost/format.hpp>
#include <boost/optional.hpp>

#include "Tags.h"
#include "WalkerPrint.h"
#include "WalkerDriver.h"
#include "Walker/InputDeck/Parser.h"
#include "Walker/CmdLine/CmdLine.h"
#include "Walker/InputDeck/InputDeck.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wconversion"
#endif

#include "distributor.decl.h"

#if defined(__clang__) || defined(__GNUC__)
  #pragma GCC diagnostic pop
#endif

namespace walker {

extern ctr::InputDeck g_inputdeck;
extern CProxy_Distributor g_DistributorProxy;

} // walker::

using walker::WalkerDriver;

WalkerDriver::WalkerDriver( const WalkerPrint& print,
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
