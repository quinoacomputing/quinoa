// *****************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     TestU01 random number generator test suite
  \details   This file declares the TestU01 random number generator test suite,
    which facilitates subjecting any supported random number generator to a
    battery of statistical tests interfaced to the TestU01 library.
*/
// *****************************************************************************
#ifndef TestU01Suite_h
#define TestU01Suite_h

#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <iosfwd>
#include <cstddef>
#include <type_traits>

#include "Tags.h"
#include "Types.h"
#include "Exception.h"
#include "StatTest.h"
#include "RNGTestPrint.h"
#include "RNGTest/InputDeck/InputDeck.h"
#include "RNGTest/Options/Battery.h"
#include "NoWarning/testu01suite.decl.h"

extern "C" {
  #include <swrite.h>
  #include <gdef.h>
}

namespace rngtest {

extern ctr::InputDeck g_inputdeck;

//! \brief TestU01 random number generator test suite used polymorphically with
//!   Battery
//! \details This class is a Charm++ chare and does the asynchronous
//!   distribution, run, and evaluation of all statistical tests within a
//!   TestU01 battery.
//! \author J. Bakosi
class TestU01Suite : public CBase_TestU01Suite {

  public:
    using Proxy = CProxy_TestU01Suite;

    //! \brief Configure all instances of the TestU01 library running on a node
    //! \details The Charm++ runtime system executes the nodeInit routine below
    //!   exactly once on every logical node early on in the Charm++ init
    //!   sequence. Must be static as it is called without and object. See also:
    //!   Section "Initializations at Program Startup" at in the Charm++ manual
    //!   http://charm.cs.illinois.edu/manuals/html/charm++/manual.html.
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
    //! \return Test suite name
    template< class Suite >
    std::string addTests() {
      const auto rngs = g_inputdeck.get< tag::selected, tag::rng >();
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

    //! Information bundle for a failed test
    struct Failed {
      std::string test;         //!< Test name
      std::string rng;          //!< RNG tested
      std::string pval;         //!< Resulting p-value
      //! Constructor with moving test name and pval string, copy RNG name
      Failed( std::string t, std::string r, std::string p ) :
        test( std::move(t) ), rng( r ), pval( std::move(p) ) {}
    };
    std::vector< Failed > m_failed;    //!< Details of failed tests
};

} // rngtest::

#endif // TestU01Suite_h
