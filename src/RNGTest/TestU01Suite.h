//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 10:25:39 AM MDT
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

    //! Constructor
    explicit TestU01Suite( ctr::BatteryType suite ) {
      // Add statistical tests to suite
      if ( suite == ctr::BatteryType::SMALLCRUSH ) addTests< SmallCrush >();
      else if ( suite == ctr::BatteryType::CRUSH ) addTests< Crush >();
      else if ( suite == ctr::BatteryType::BIGCRUSH ) addTests< BigCrush >();
      else Throw( tk::ExceptType::FATAL,
                  "Non-TestU01 RNG test suite passed to TestU01Suite" );

      // Call test constructors and fire up collection of number of results and
      // test names as futures collected in battery()
      std::size_t numtests = ntest();
      std::size_t n = 0;
      for (const auto& t : m_tests) {
        const auto& test = t();
        if (n++ < numtests) {      // only collect from one test per RNGs tested
          f_nstat.emplace_back( CkCreateFuture() );
          test.npval( f_nstat.back() );
          f_name.emplace_back( CkCreateFuture() );
          test.name( f_name.back() );
        }
      }
    }

    //! Run battery of RNG tests
    void run() {
      print();
      //for (const auto& t : m_tests) t().run(0);
      mainProxy.finalize();
    }

    //! Evaluate a statistical test
    void evaluate( std::size_t id ) const {
      std::cout << "eval: " << id << std::endl;
      mainProxy.finalize();
    }

  private:
    //! Print information on test suite
    void print() const {
      RNGTestPrint print;
      battery( print );
      print.section( "RNGs tested" );
      #ifdef HAS_MKL
      print.MKLParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                       g_inputdeck.get< tag::param, tk::tag::rngmkl >() );
      #endif
      print.RNGSSEParams( g_inputdeck.get< tag::selected, tk::tag::rng >(),
                          g_inputdeck.get< tag::param, tk::tag::rngsse >() );
      print.raw('\n');
    }

    //! Add all statistical tests to suite
    template< class Suite >
    void addTests() {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      ErrChk( !rngs.empty(), tk::ExceptType::FATAL, "No RNGs selected" );
      Suite suite;
      for (const auto& s : rngs) suite.addTests( m_tests, s, thisProxy );
    }

    //! Return number of statistical tests
    std::size_t ntest() const {
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      return m_tests.size() / rngs.size();
    }

    //! Print information on battery
    void battery( RNGTestPrint& print ) const {
      std::size_t nstat = 0;
      for (const auto& n : f_nstat)
        nstat += tk::waitfor< tk::Msg<std::size_t> >( n );
      const auto rngs = g_inputdeck.get< tag::selected, tk::tag::rng >();
      print.battery( m_tests.size()/rngs.size(), nstat );
      for (const auto& n : f_name)
        print.names( tk::waitfor< tk::StringsMsg >( n ) );
    }

    std::vector< std::function< StatTest() > > m_tests; //!< Statistical tests
    std::vector< CkFuture > f_nstat;   //!< Futures collecting number of results
    std::vector< CkFuture > f_name;    //!< Futures collecting test names
};

} // rngtest::

#endif // TestU01Suite_h
