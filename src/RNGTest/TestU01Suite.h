//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Fri 25 Jul 2014 09:40:44 AM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
#include <testu01suite.decl.h>

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
    explicit TestU01Suite( ctr::BatteryType suite );

    //! Collect number of p-values from a test
    void npval( std::size_t n );

    //! Collect test name(s) from a test
    void names( std::vector< std::string > n );

    //! Evaluate a statistical test
    void evaluate( std::vector< std::vector< std::string > > status );

    //! Collect test run time from a test
    void time( std::pair< std::string, tk::real > t );

 private:
    //! Add all statistical tests to suite, return suite name
    template< class Suite >
    std::string addTests() {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      ErrChk( !rngs.empty(), "No RNGs selected" );
      Suite suite;
      for (const auto& r : rngs) suite.addTests( m_ctrs, r, thisProxy );
      return suite.name();
    }

    //! Return number of statistical tests
    std::size_t ntest() const;

    //! Output final assessment
    void assess();

    RNGTestPrint m_print;              //!< Pretty printer
    std::vector< std::function< StatTest() > > m_ctrs; //! Tests constructors
    std::vector< StatTest > m_tests;   //!< Constructed statistical tests
    std::string m_name;                //!< Test suite name
    std::size_t m_npval;               //!< Number of results from all tests
    std::size_t m_ncomplete;           //!< Number of completed tests
    std::size_t m_ntest;               //!< Number of tests info received from
    std::map< std::string, std::size_t > m_nfail; //! Number of failed tests/RNG
    std::map< std::string, tk::real > m_time;     //!< Measured time/RNG

    struct Failed {
      std::string test;
      std::string rng;
      std::string pval;
      // Construct with moving test name and pval string, copy RNG name
      Failed( std::string t, std::string r, std::string p ) :
        test( std::move(t) ), rng( r ), pval( std::move(p) ) {}
    };
    std::vector< Failed > m_failed;    //!< Details of failed tests
};

} // rngtest::

#endif // TestU01Suite_h
