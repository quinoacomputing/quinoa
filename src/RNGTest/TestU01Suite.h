//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Sat 28 Jun 2014 09:57:46 PM MDT
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
#include <Msg.h>
#include <testu01suite.decl.h>
#include <rngtest.decl.h>       // only for mainProxy, will go away

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
    explicit TestU01Suite( ctr::BatteryType suite ) : m_npval(0), m_ncomplete(0)
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

      std::vector< StatTest > tests;
      // Construct tests and start collecting number of results and test names
      std::size_t numtests = ntest();
      std::size_t n = 0;
      for (const auto& t : m_tests) {
        tests.push_back( t() );         // construct test and store handle
        if (n++ < numtests) {           // only collect one test per RNGs tested
          //// Fire up collection of number of results
          f_npval.emplace_back( CkCreateFuture() );
          tests.back().npval( f_npval.back() );
          // Fire up collection of test names
          f_name.emplace_back( CkCreateFuture() );
          tests.back().name( f_name.back() );
        }
      }

      // Print out info on test suite
      print();

      // Run battery of RNG tests
      for (const auto& t : tests) t.run();

      // Initialize space for counting the number failed tests per RNG. Note
      // that we could use tk::ctr::RNGType as the map-key here instead of the
      // RNG name (RNGType would be an enum vs the name which is a string), but
      // then the RNGType would have to be part of the status (instead of the
      // RNG name) returning from a completed test by TestU01Props::run(). That
      // would require a return type that is less generic than the current
      // vector of vector of strings. Since the RNG name is already part of the
      // status, we just match the RNG name for counting failed tests per RNG
      // and keep the status more generic. See also the discussion on the return
      // type in TestU01Props::run().
      tk::Option< tk::ctr::RNG > rng;
      for (const auto& r : g_inputdeck.get< tag::selected, tk::tag::rng >() )
        m_nfail[ rng.name(r) ] = 0;
    }

    //! Evaluate a statistical test
    void evaluate( const std::vector< std::vector< std::string > >& status ) {
      m_print.test( ++m_ncomplete, m_tests.size(), m_nfail, status );
      if ( m_ncomplete == m_tests.size() ) mainProxy.finalize();
    }

  private:
    //! Print information on test suite
    void print() {
      // Wait for and add up total number of p-values expected
      for (const auto& n : f_npval)
        m_npval += tk::waitfor< tk::Msg<std::size_t> >( n );
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      m_print.battery( m_tests.size()/rngs.size(), m_npval );
      // Wait for and print the test names
      for (const auto& n : f_name)
        m_print.names( tk::waitfor< tk::StringsMsg >( n ) );

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
    }

    //! Add all statistical tests to suite, return suite name
    template< class Suite >
    std::string addTests() {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      ErrChk( !rngs.empty(), tk::ExceptType::FATAL, "No RNGs selected" );
      Suite suite;
      for (const auto& s : rngs) suite.addTests( m_tests, s, thisProxy );
      return suite.name();
    }

    //! Return number of statistical tests
    std::size_t ntest() const {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      return m_tests.size() / rngs.size();
    }

    RNGTestPrint m_print;               //!< Pretty printer
    std::vector< std::function< StatTest() > > m_tests; //!< Statistical tests
    std::vector< CkFuture > f_npval;   //!< Futures collecting number of results
    std::vector< CkFuture > f_name;    //!< Futures collecting test names
    std::string m_name;                //!< Test suite name
    std::size_t m_npval;               //!< Number of results from all tests
    std::size_t m_ncomplete;           //!< Number of completed tests
    //! Number of failed tests per RNG
    std::map< std::string, std::size_t > m_nfail;
};

} // rngtest::

#endif // TestU01Suite_h
