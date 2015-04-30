//******************************************************************************
/*!
  \file      src/Main/RNGTestDriver.C
  \author    J. Bakosi
  \date      Fri 17 Apr 2015 11:50:23 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Random number generator test suite driver
  \details   Random number generator test suite driver.
*/
//******************************************************************************

#include <Factory.h>
#include <RNGTestDriver.h>
#include <RNGTest/InputDeck/Parser.h>
#include <TestU01Suite.h>

#ifdef HAS_MKL
#include <MKLRNG.h>
#endif

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

} // rngtest::

using rngtest::RNGTestDriver;

RNGTestDriver::RNGTestDriver( const RNGTestPrint& print,
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
  // All global-scope data to be migrated to all PEs initialized here

  // Parse input deck into g_inputdeck, transfer cmdline (no longer needed)
  m_print.item( "Control file", cmdline.get< tag::io, tag::control >() );  
  InputDeckParser inputdeckParser( m_print, cmdline, g_inputdeck );
  m_print.item( "Parsed control file", "success" );  
  m_print.endpart();
}

void
RNGTestDriver::execute() const
//******************************************************************************
//  Run battery
//! \author J. Bakosi
//******************************************************************************
{
  m_print.part( "Factory" );

  // Register batteries
  BatteryFactory bf;
  using ctr::BatteryType;
  // Note that TestU01Suite constructors take the BatteryType (enum class)
  // value as their argument, which happens to be the same as the key in the
  // factory - hence the double-specification of the battery type below.
  tk::recordCharmModel< Battery, TestU01Suite >
                      ( bf, BatteryType::SMALLCRUSH, BatteryType::SMALLCRUSH );
  tk::recordCharmModel< Battery, TestU01Suite >
                      ( bf, BatteryType::CRUSH, BatteryType::CRUSH );
  tk::recordCharmModel< Battery, TestU01Suite >
                      ( bf, BatteryType::BIGCRUSH, BatteryType::BIGCRUSH );
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
