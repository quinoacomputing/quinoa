//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.C
  \author    J. Bakosi
  \date      Fri 01 Aug 2014 11:42:36 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     TestU01 suite
  \details   TestU01 suite
*/
//******************************************************************************

#include <TestU01Suite.h>
#include <TestStack.h>
#include <RNGTest/Options/Battery.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>
#include <rngtest.decl.h>

extern CProxy_Main mainProxy;

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::TestU01Suite;

TestU01Suite::TestU01Suite( ctr::BatteryType suite ) :
  m_print( rngtest::g_inputdeck.get< rngtest::tag::cmd, tk::tag::verbose >() ?
           std::cout : std::clog ),
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
  else Throw( "Non-TestU01 RNG test suite passed to TestU01Suite" );

  // Construct all tests and store handles
  for (const auto& t : m_ctrs) m_tests.emplace_back( t() );

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
    const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
    std::stringstream ss;
    ss << "RNGs tested (" << rngs.size() << ")";
    m_print.section( ss.str() );
    #ifdef HAS_MKL
    m_print.MKLParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                       g_inputdeck.get< tag::param, tk::tag::rngmkl >() );
    #endif
    m_print.RNGSSEParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                          g_inputdeck.get< tag::param, tk::tag::rngsse >() );
    m_print.endpart();
    m_print.part( m_name );
    m_print.statshead( "Statistics computed",
                       m_npval*rngs.size(),
                       m_ctrs.size() );

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
  m_print.test( ++m_ncomplete, m_ctrs.size(), m_nfail, status );

  // Store information on failed test for final assessment
  for (std::size_t p=0; p<status[1].size(); ++p)
    if (status[1][p].size() > 4)
      m_failed.emplace_back( status[0][p], status[2][0], status[1][p] );

  if ( m_ncomplete == m_ctrs.size() ) {
    // Collect measured test run times
    m_ncomplete = 0;
    for (const auto& t : m_tests) t.time();
  }
}

void
TestU01Suite::time( std::pair< std::string, tk::real > t )
//******************************************************************************
// Collect test times measured in seconds from a statistical test
//! \author  J. Bakosi
//******************************************************************************
{
  m_time[ t.first ] += t.second;

  if ( ++m_ncomplete == m_ctrs.size() ) assess();
}

void
TestU01Suite::assess()
//******************************************************************************
// Output final assessment
//! \author  J. Bakosi
//******************************************************************************
{
  // Output summary of failed tests for all RNGs tested
  if ( !m_failed.empty() ) {
    const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
    m_print.failed( "Failed statistics", m_npval*rngs.size(), m_failed );
  } else m_print.note< tk::QUIET >( "All tests passed" );

  // Cost and quality assessment only for more than one RNG
  if (m_time.size() > 1) {
    // Output measured times per RNG in order of computational cost
    m_print.cost( "Generator cost",
                  "Measured times in seconds in increasing order (low is good)",
                  m_time );
    // Output number of failed tests per RNG in order of decreasing quality
    m_print.rank( "Generator quality",
                  "Number of failed tests in increasing order (low is good)",
                  m_nfail );
  }

  // Quit
  mainProxy.finalize();
}

std::size_t
TestU01Suite::ntest() const
//******************************************************************************
// Return the number of statistical tests per each RNG tested
//! \author  J. Bakosi
//******************************************************************************
{
  const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
  return m_ctrs.size() / rngs.size();
}

#include <testu01suite.def.h>
