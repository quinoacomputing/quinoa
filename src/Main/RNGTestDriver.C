// *****************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Random number generator test suite driver
  \details   Random number generator test suite driver.
*/
// *****************************************************************************

#include <string>
#include <utility>
#include <iterator>

#include "Tags.h"
#include "Exception.h"
#include "Factory.h"
#include "Battery.h"
#include "TestU01Suite.h"
#include "RNGTestPrint.h"
#include "RNGTestDriver.h"
#include "RNGTest/InputDeck/InputDeck.h"
#include "RNGTest/InputDeck/Parser.h"

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

} // rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver( const RNGTestPrint& print,
                              const ctr::CmdLine& cmdline ) :
  m_print( print )
// *****************************************************************************
//  Constructor
//! \param[in] print Pretty printer
//! \param[in] cmdline Command line object storing data parsed from the command
//!   line arguments
//! \author J. Bakosi
// *****************************************************************************
{
  // All global-scope data to be migrated to all PEs initialized here

  // Parse input deck into g_inputdeck, transfer cmdline (no longer needed)
  m_print.item( "Control file", cmdline.get< tag::io, tag::control >() );  
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );
  m_print.item( "Parsed control file", "success" );  
  m_print.endpart();
}

void
RNGTestDriver::execute() const
// *****************************************************************************
//  Run battery
//! \author J. Bakosi
// *****************************************************************************
{
  m_print.part( "Factory" );

  // Register batteries
  BatteryFactory bf;
  using ctr::BatteryType;
  // Note that TestU01Suite constructors take the BatteryType (enum class)
  // value as their argument, which happens to be the same as the key in the
  // factory - hence the double-specification of the battery type below.
  // Record all into a factory passing the last 0 means instantiate on PE 0.
  tk::recordCharmModel< Battery, TestU01Suite >
                    ( bf, BatteryType::SMALLCRUSH, BatteryType::SMALLCRUSH, 0 );
  tk::recordCharmModel< Battery, TestU01Suite >
                    ( bf, BatteryType::CRUSH, BatteryType::CRUSH, 0 );
  tk::recordCharmModel< Battery, TestU01Suite >
                    ( bf, BatteryType::BIGCRUSH, BatteryType::BIGCRUSH, 0 );
  m_print.list< ctr::Battery >( "Registered batteries", bf );
  m_print.endpart();

  m_print.part( "Problem" );
  if ( !g_inputdeck.get< tag::title >().empty() )
    m_print.title( g_inputdeck.get< tag::title >() );

  // Instantiate and run battery
  const auto s = bf.find( g_inputdeck.get< tag::selected, tag::battery >() );
  if (s != end(bf)) {
    s->second();
  } else Throw( "Battery not found in factory" );
}
