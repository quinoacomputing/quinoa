//******************************************************************************
/*!
  \file      src/RNGTest/TestU01Suite.h
  \author    J. Bakosi
  \date      Tue 03 Jun 2014 09:09:26 PM MDT
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

#include <array>

#include <Base.h>
#include <Options/RNG.h>
#include <PUPUtil.h>
#include <Battery.h>
#include <TestU01Util.h>
#include <testu01suite.decl.h>

namespace rngtest {

extern const std::size_t g_maxRNGs;
extern int g_tid;

//! TestU01 random number generator test suite
//! Suite is a policy specific to a test suite, see use cases below
template< class Suite >
class TestU01Suite : public Battery {

  public:
    PUPable_decl_template( TestU01Suite );

    //! Constructor
    explicit TestU01Suite( const Base& base ) :
      m_numRNGs( 1 ),
      m_rng( base.control.get< tag::selected, tk::tag::rng >() ),
      m_npval(0),
      m_ncomplete( 0 )
    {
      Suite suite;
      m_name = suite.policy();
      //assignTests();
      PUPable_reg( TestU01Suite );
    }

    //! Destructor
    ~TestU01Suite() override = default;

    //! Migrator
    TestU01Suite( CkMigrateMessage* m = 0 ) {}

//     // Add a suite of tests for each RNG tested
//     void assignTests() {
//       for (std::size_t sel=0; sel<m_rng.size(); ++sel)
//         suite.addTests( sel, m_rng, m_tests );
//     }

    //! Pack/Unpack statistical tests polymorfic pointers
    void pup_tests( PUP::er& p ) {
      auto size = tk::pup_container_size( p, m_tests );
      if ( p.isUnpacking() ) m_tests.resize( size );
      //assignTests();
    }

    //! Pack/Unpack
    void pup( PUP::er& p ) override {
      //pup_vector( p, m_rng );
      p | m_npval;
      p | m_ncomplete;
      p | m_name;
      pup_tests( p );
    }

    //! Run battery of RNG tests
    void run() override {
      std::cout << "TestU01Suite::run()" << std::endl;
//  const RNGTestPrint& pr = *g_base.print;
//  const tk::Timer& timer = *g_base.timer;

//   pr.part( m_name );
//   pr.statshead( "Statistics computed" );

      swrite_Basic = FALSE;         // want no screen putput from TestU01

      std::size_t nt = m_tests.size();
      //std::array< std::size_t, g_maxRNGs > nfail{{ 0 }};

      #ifdef _OPENMP
      #pragma omp parallel
      #endif
      {
        #ifdef _OPENMP
        //g_tid = omp_get_thread_num();
        #else
        //g_tid = 0;
        #endif

        #ifdef _OPENMP
        #pragma omp for schedule(dynamic)
        #endif
        for ( std::size_t i=0; i<nt; ++i ) {

          // Get reference to ith test
          std::unique_ptr< StatTest >& test = m_tests[i];

          // Get RNG id being tested
          auto id = test->id();

          // Get RNG properties
          //auto& props = m_rngprops[id];

//       // Start timer for RNG
//       timer.start( props.timer[g_tid] );

      // Run test
      //test->run();
//      CProxy_StatTestRunner::ckNew( *test, i );

//       // Query and accumulate elapsed time for RNG
//       tk::real time = timer.query( props.timer[g_tid] );
//       #ifdef _OPENMP
//       #pragma omp atomic
//       #endif
//       props.time[g_tid] += time;
// 
//       // Evaluate test
//       auto npval = test->nstat();
//       for (std::size_t p=0; p<npval; ++p) {
// 
// //         // Increase number tests completed
// //         #ifdef _OPENMP
// //         #pragma omp atomic
// //         #endif
// //         ++ncomplete;
// 
//         // Increase number of failed tests for RNG
//         if ( test->fail(p) ) {
//           #ifdef _OPENMP
//           #pragma omp atomic
//           #endif
//           ++nfail[id];
//         }
// 
//         // Output one-liner
//         pr.test< StatTest, TestContainer >
//                ( m_ncomplete, nfail[id], m_npval, test, p );
//       }
        }
      }

//   // Count up number of total failed tests (for all RNGs tested)
//   auto tfail = failed();
// 
//   // Output summary of failed tests (for all RNGs tested)
//   if (tfail) {
//     pr.failed< StatTest >("Failed statistics", m_npval, tfail, m_tests);
//   } else {
//     pr.note("All tests passed");
//   }
// 
//   // Compute sum of measured times spent by all threads per RNG
//   std::vector< std::pair< tk::real, std::string > > rngtimes;   // times & names
//   std::vector< std::pair< std::size_t, std::string > > rngnfail;// nfail & names
//   tk::ctr::RNG rng;
//   std::size_t i=0;
//   for (const auto& r : m_rngprops) {
//     rngtimes.push_back( { 0.0, rng.name(r.id) } );
//     rngnfail.push_back( { nfail[i], rng.name(r.id) } );
//     for (const auto& t : r.time) {
//       rngtimes.back().first += t;
//     }
//     if ( std::fabs(rngtimes.back().first) <     // remove if not tested
//          std::numeric_limits<tk::real>::epsilon() ) {
//       rngtimes.pop_back();
//       rngnfail.pop_back();
//     }
//     ++i;
//   }
// 
//   // Output measured times per RNG in order of computational cost
//   pr.cost( "Generator cost",
//            "Measured times in seconds in increasing order (low is good)",
//            rngtimes );
// 
//   // Output number of failed tests per RNG in order of decreasing quality
//   pr.rank( "Generator quality",
//            "Number of failed tests in increasing order (low is good)",
//            rngnfail );
// 
//   pr.endpart();
    }

    //! Print list of registered statistical tests
    void print( const RNGTestPrint& print ) const override {
     // Output test names (only for the first RNG tested, the rest are repeats)
     print.names< StatTest >( m_tests, ntest() );
    }

    //! Return number of statistical tests in battery
    std::size_t ntest() const override {
      return m_numRNGs ? m_tests.size() / m_numRNGs : 0;
    }

    //! Return number of statistics produced by battery
    std::size_t nstat() const override {
      return m_numRNGs ? m_npval / m_numRNGs : 0;
    }

    void evaluate( std::size_t id ) override {
      std::cout << "TestU01Suite::evaluate( " << id << " )" << std::endl;
      std::cout << "m_name.size(): " << m_name.size() << std::endl;
      //std::cout << "ntest: " << ntest() << std::endl;
      if ( ++m_ncomplete == 1/*ntest()*/ ) {
        CkPrintf("All tests evaluated.\n");
        CkExit();
      } else {
        CkPrintf("Waiting for more tests.\n");
      }
    }

  private:
    //! Don't permit copy constructor
    TestU01Suite(const TestU01Suite&) = delete;
    //! Don't permit copy assigment
    TestU01Suite& operator=(const TestU01Suite&) = delete;
    //! Don't permit move constructor
    TestU01Suite(TestU01Suite&&) = delete;
    //! Don't permit move assigment
    TestU01Suite& operator=(TestU01Suite&&) = delete;

    //! Count up total number of statistics expected
    void total() {
      m_npval = 0;
      for (auto& t : m_tests) m_npval += t->nstat();
    }

    //! Count up number of failed tests
    std::size_t failed() {
      std::size_t fail = 0;
      for (auto& t : m_tests) fail += t->nfail();
      return fail;
    }

    std::size_t m_numRNGs;
    std::vector< ctr::RNGType > m_rng;  //!< Selected RNGs (from control)
    std::size_t m_npval;                //!< Total number of stats
    std::size_t m_ncomplete;            //!< Count number of completed tests
    std::string m_name;                 //!< Test suite name
    std::vector< std::unique_ptr< StatTest > > m_tests; //!< Statistical tests
};

} // rngtest::

#endif // TestU01Suite_h
