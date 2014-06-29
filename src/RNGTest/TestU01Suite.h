//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Sun 29 Jun 2014 05:07:33 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     TestU01 random number generator test suite
  \details   TestU01 random number generator test suite
*/
//******************************************************************************
#ifndef TestU01Suite_h
#define TestU01Suite_h

extern "C" {
  #include <swrite.h>
}

#include <StatTest.h>
#include <RNGTestPrint.h>
#include <RNGTest/Options/Battery.h>
#include <SmallCrush.h>
#include <Crush.h>
#include <BigCrush.h>
#include <testu01suite.decl.h>
#include <rngtest.decl.h>

extern CProxy_Main mainProxy;

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

//! TestU01 random number generator test suite used polymorphically with Battery
class TestU01Suite : public CBase_TestU01Suite {

  public:
    using Proxy = CProxy_TestU01Suite;

    //! The Charm++ runtime system executes the nodeInit routine below exactly
    //! once on every logical node early on in the Charm++ init sequence. Must
    //! be static as it is called without and object. See also: Section
    //! "Initializations at Program Startup" at http://charm.cs.illinois.edu/
    //! manuals/html/charm++/manual.html.
    static void nodeInit() {
      swrite_Basic = FALSE;     // want no screen output from TestU01 library
    }

    //! Constructor
    explicit TestU01Suite( ctr::BatteryType suite ) :
      m_npval(0), m_ncomplete(0), m_ntest(0)
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

    //! Collect number of p-values from a test
    void npval( std::size_t n ) {
      m_npval += n;
      if ( ++m_ntest == ntest() ) {
        m_print.battery( ntest(), m_npval );
        // Collect test names from all tests (one per RNG)
        m_ntest = 0;
        for (std::size_t i=0; i<ntest(); ++i) m_tests[i].names();
      }
    }

    //! Collect test name(s) from a test
    void names( std::vector< std::string > n ) {
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

    //! Evaluate a statistical test
    void evaluate( std::vector< std::vector< std::string > > status ) {
      m_print.test( ++m_ncomplete, m_testFactory.size(), m_nfail, status );
      if ( m_ncomplete == m_testFactory.size() ) mainProxy.finalize();
    }

 private:
    //! Add all statistical tests to suite, return suite name
    template< class Suite >
    std::string addTests() {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      ErrChk( !rngs.empty(), tk::ExceptType::FATAL, "No RNGs selected" );
      Suite suite;
      for (const auto& s : rngs) suite.addTests( m_testFactory, s, thisProxy );
      return suite.name();
    }

    //! Return number of statistical tests
    std::size_t ntest() const {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      return m_testFactory.size() / rngs.size();
    }

    RNGTestPrint m_print;               //!< Pretty printer
    //! Bound statistical tests constructors
    std::vector< std::function< StatTest() > > m_testFactory;
    std::vector< StatTest > m_tests;   //!< Constructed statistical tests
    std::string m_name;                //!< Test suite name
    std::size_t m_npval;               //!< Number of results from all tests
    std::size_t m_ncomplete;           //!< Number of completed tests
    std::size_t m_ntest;               //!< Number of tests info received from
    //! Number of failed tests per RNG
    std::map< std::string, std::size_t > m_nfail;
};

} // rngtest::

#endif // TestU01Suite_h
