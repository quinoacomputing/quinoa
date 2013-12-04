//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Wed 04 Dec 2013 03:43:00 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <iostream>

extern "C" {
  #include <unif01.h>
  #include <bbattery.h>
  #include <swrite.h>
}

#include <SmallCrush.h>
#include <TestU01.h>

namespace rngtest {

extern int g_tid;               //!< Global thread id (in TestU01Suite.C)
#ifdef _OPENMP
#pragma omp threadprivate(g_tid)
#endif

} // rngtest::

using rngtest::SmallCrush;

SmallCrush::SmallCrush(const Base& base) : TestU01Suite(base)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Get vector of selected RNGs
  std::vector< ctr::RNGType > rngs =
    m_base.control.get< ctr::selected, ctr::rng >();

  // Instantiate all selected RNGs and add statistical tests for each
  for (const auto& r : rngs) {
    if (r != ctr::RNGType::NO_RNG) {
      // ICC: replace linear search with map.find()
      // const auto& rng = m_rng.find(r);   // Find RNG in registry
      // Assert( rng != m_rng.end(), tk::ExceptType::FATAL, "RNG not found" );
      // m_testRNGs.push_back( rng->first );// Build vector of RNG ids to test
      // addTests( r, rng->second );        // Add tests to battery for this RNG
      for (size_t i=0; i<m_rng1.size(); ++i) {
        if (m_rng1[i] == r) {
          m_testRNGs.push_back( r );        // Build vector of RNG ids to test
          addTests( r, m_rng[i] );          // Add tests to battery for this RNG
        }
      }
    }
  }

  // Count up total number of statistics expected
  m_npval = 0;
  for (const auto& t : m_pvals) {
    m_npval += t.size();
  }
}

void
SmallCrush::addTests( const quinoa::ctr::RNGType& rng, const Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
  m_pvals.push_back(
    add< TestU01< sres_Poisson,
                  sres_CreatePoisson,
                  sres_DeletePoisson,
                  BirthdaySpacings > >
       ( m_tests, gen, rng, StatTest::Names({"Birthday Spacings"}) ) );
  
  m_pvals.push_back(
    add< TestU01< sknuth_Res2,
                  sknuth_CreateRes2,
                  sknuth_DeleteRes2,
                  Collision > >
       ( m_tests, gen, rng, StatTest::Names({"Collision"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  Gap > >
       ( m_tests, gen, rng, StatTest::Names({"Gap"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  SimpPoker > >
       ( m_tests, gen, rng, StatTest::Names({"Simplified Poker"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  CouponCollector > >
       ( m_tests, gen, rng, StatTest::Names({"Coupon Collector"}) ) );

  m_pvals.push_back(
    add< TestU01< sknuth_Res1,
                  sknuth_CreateRes1,
                  sknuth_DeleteRes1,
                  MaxOft > >
       ( m_tests, gen, rng,
         StatTest::Names({"Maximum-of-t",
                          "Maximum-of-t Anderson-Darling"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  WeightDistrib > >
       ( m_tests, gen, rng,
         StatTest::Names({"Weight Distribution"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  MatrixRank > >
       ( m_tests, gen, rng, StatTest::Names({"Matrix Rank"}) ) );

  m_pvals.push_back(
    add< TestU01< sstring_Res,
                  sstring_CreateRes,
                  sstring_DeleteRes,
                  HammingIndep > >
       ( m_tests, gen, rng,
         StatTest::Names({"Hamming Independence"}) ) );

  m_pvals.push_back(
    add< TestU01< swalk_Res,
                  swalk_CreateRes,
                  swalk_DeleteRes,
                  RandomWalk1 > >
       ( m_tests, gen, rng,
         StatTest::Names({"Random Walk 1 Stat H",
                          "Random Walk 1 Stat M",
                          "Random Walk 1 Stat J",
                          "Random Walk 1 Stat R",
                          "Random Walk 1 Stat C"}) ) );
}

void
SmallCrush::run()
//******************************************************************************
//  Run SmallCrush battery
//! \author  J. Bakosi
//******************************************************************************
{
  const RNGTestPrint& print = m_base.print;

  print.part("Run SmallCrush");
  print.section("Statistical tests finished");

  swrite_Basic = FALSE;         // Want screen no putput from TestU01

//   g_tid = 0;
//   bbattery_SmallCrush( m_rng[2].get() );

  using Pval = StatTest::Pvals::value_type;
  using Psize = StatTest::Pvals::size_type;
  using Tsize = TestContainer::size_type;

  Tsize ntest = m_tests.size();
  Psize n = 0;
  #ifdef _OPENMP
  #pragma omp parallel
  #endif
  {
    #ifdef _OPENMP
    g_tid = omp_get_thread_num();
    #else
    g_tid = 0;
    #endif

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (Tsize i=0; i<ntest; ++i) {

      // Run test
      StatTest::Pvals pvals = m_tests[i]->run();

      // Evaluate test
      Psize npval = pvals.size();
      for (Psize p=0; p<npval; ++p) {

        // Increase number tests completed
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        ++n;

        // Output one-liner
        m_base.print.test< StatTest, TestContainer >
                         ( n, m_npval, m_tests[i], p, pvals[p] );

        // Save suspect p-value
        if ((pvals[p] <= gofw_Suspectp) || (pvals[p] >= 1.0 - gofw_Suspectp)) {
          m_pvals[i][p] = pvals[p];
        }
      }
    }
  }

  // Count up number of failed tests
  StatTest::Pvals::size_type failed = 0;
  for (const auto& t : m_pvals) {
    for (const auto& p : t) {
      if (fabs(p+1.0) > std::numeric_limits< Pval >::epsilon()) {
        ++failed;
      }
    }
  }

  // Output failed tests
  if (failed) {
    m_base.print.failed< StatTest >
                       ( "Failed tests", m_npval, failed, m_pvals, m_tests );
  } else {
    m_base.print.note("All tests passed");
  }

  m_base.print.endpart();
}

void
SmallCrush::print() const
//******************************************************************************
//  Print list of registered statistical tests
//! \author  J. Bakosi
//******************************************************************************
{
  // Output test names (only for the first RNG tested, the rest are repeats)
  m_base.print.names< StatTest >( m_tests, m_tests.size()/m_testRNGs.size() );
}
