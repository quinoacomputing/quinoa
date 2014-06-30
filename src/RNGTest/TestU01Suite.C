//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.C
  \author    J. Bakosi
  \date      Mon 30 Jun 2014 07:00:26 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

#include <TestU01Suite.h>
#include <RNGTest/Options/Battery.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>
#include <rngtest.decl.h>

extern CProxy_Main mainProxy;

using rngtest::TestU01Suite;

TestU01Suite::TestU01Suite( ctr::BatteryType suite ) :
  m_npval(0), m_ncomplete(0), m_ntest(0)
//******************************************************************************
// Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Add statistical tests to suite
  if ( suite == ctr::BatteryType::SMALLCRUSH )
    m_name = addTests< SmallCrush >();
  else if ( suite == ctr::BatteryType::CRUSH )
    m_name = addTests< Crush >();
  else if ( suite == ctr::BatteryType::BIGCRUSH )
    m_name = addTests< BigCrush >();
  else Throw( tk::ExceptType::FATAL,
              "Non-TestU01 RNG test suite passed to TestU01Suite" );

  // Construct all tests and store handles
  for (const auto& t : m_testFactory) m_tests.emplace_back( t() );

  // Collect number of results from all tests (one per RNG)
  for (std::size_t i=0; i<ntest(); ++i) m_tests[i].npval();
}

void
TestU01Suite::npval( std::size_t n )
//******************************************************************************
// Collect number of p-values from a statistical test
//! \author  J. Bakosi
//******************************************************************************
{
  m_npval += n;
  if ( ++m_ntest == ntest() ) {
    m_print.battery( ntest(), m_npval );
    // Collect test names from all tests (one per RNG)
    m_ntest = 0;
    for (std::size_t i=0; i<ntest(); ++i) m_tests[i].names();
  }
}

void
TestU01Suite::names( std::vector< std::string > n )
//******************************************************************************
// Collect test names from a statistical test
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.names( n );
  if ( ++m_ntest == ntest() ) {
    m_print.section( "RNGs tested" );
    #ifdef HAS_MKL
    m_print.MKLParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                       g_inputdeck.get< tag::param, tk::tag::rngmkl >() );
    #endif
    m_print.RNGSSEParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                          g_inputdeck.get< tag::param, tk::tag::rngsse >() );
    m_print.raw('\n');
    m_print.part( m_name );
    m_print.statshead( "Statistics computed" );

    // Run battery of RNG tests
    for (const auto& t : m_tests) t.run();

    // Initialize space for counting the number failed tests per RNG. Note
    // that we could use tk::ctr::RNGType as the map-key here instead of the
    // RNG name (RNGType would be an enum vs the name which is a string),
    // but then the RNGType would have to be part of the status (instead of
    // the RNG name) returning from a completed test by TestU01Props::run().
    // That would require a return type that is less generic than the
    // current vector of vector of strings. Since the RNG name is already
    // part of the status, we just match the RNG name for counting failed
    // tests per RNG and keep the status more generic. See also the
    // discussion on the return type in TestU01Props::run().
    tk::Option< tk::ctr::RNG > rng;
    for (const auto& r : g_inputdeck.get< tag::selected, tk::tag::rng >() )
      m_nfail[ rng.name(r) ] = 0;
  }
}

void
TestU01Suite::evaluate( std::vector< std::vector< std::string > > status )
//******************************************************************************
// Evaluate statistical test
//! \author  J. Bakosi
//******************************************************************************
{
  m_print.test( ++m_ncomplete, m_testFactory.size(), m_nfail, status );
  if ( m_ncomplete == m_testFactory.size() ) mainProxy.finalize();
}

std::size_t
TestU01Suite::ntest() const
//******************************************************************************
// Return the number of statistical tests per each RNG tested
//! \author  J. Bakosi
//******************************************************************************
{
  const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
  return m_testFactory.size() / rngs.size();
}

#include <testu01suite.def.h>
