//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Fri 29 Nov 2013 09:58:54 AM MST
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
#include <MKLRNGWrappers.h>
#include <TestU01.h>

namespace rngtest {

extern std::unique_ptr< tk::RNG > g_rng;        //!< RNG used by TestU01
extern int g_tid;                               //!< Thread id used by TestU01

} // rngtest::

using rngtest::SmallCrush;

SmallCrush::SmallCrush(const Base& base) : TestU01Suite(base)
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  // Instantiate RNGs
  ctr::RNGType r = m_base.control.get< ctr::selected, ctr::rng >()[0];
  if (r != ctr::RNGType::NO_RNG) {
    // Save name of random number generator
    quinoa::ctr::RNG rng;
    m_rngname = rng.name( r );
    // Instantiate random number generator
    g_rng = std::unique_ptr< tk::RNG >( m_base.rng[r]() );
  }

  // Create TestU01 external generator
  char* const rngname = const_cast<char*>( m_rngname.c_str() );
  m_gen = Gen01Ptr( unif01_CreateExternGen01( rngname, MKLRNGUniform ) );

  // Add statistical tests to battery
  addTests();
}

void
SmallCrush::addTests()
//******************************************************************************
// Add statistical tests to battery, count up total number of p-values
//! \author  J. Bakosi
//******************************************************************************
{
  m_pvals.push_back(
    add< TestU01< sres_Poisson,
                  sres_CreatePoisson,
                  sres_DeletePoisson,
                  BirthdaySpacings > >
       ( m_tests, m_gen, StatTest::Names({"Birthday Spacings"}) ) );
  

  m_pvals.push_back(
    add< TestU01< sknuth_Res2,
                  sknuth_CreateRes2,
                  sknuth_DeleteRes2,
                  Collision > >
       ( m_tests, m_gen, StatTest::Names({"Collision"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  Gap > >
       ( m_tests, m_gen, StatTest::Names({"Gap"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  SimpPoker > >
       ( m_tests, m_gen, StatTest::Names({"Simplified Poker"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  CouponCollector > >
       ( m_tests, m_gen, StatTest::Names({"Coupon Collector"}) ) );

  m_pvals.push_back(
    add< TestU01< sknuth_Res1,
                  sknuth_CreateRes1,
                  sknuth_DeleteRes1,
                  MaxOft > >
       ( m_tests, m_gen,
         StatTest::Names({"Maximum-of-t",
                          "Maximum-of-t Anderson-Darling"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  WeightDistrib > >
       ( m_tests, m_gen,
         StatTest::Names({"Weight Distribution"}) ) );

  m_pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  MatrixRank > >
       ( m_tests, m_gen, StatTest::Names({"Matrix Rank"}) ) );

  m_pvals.push_back(
    add< TestU01< sstring_Res,
                  sstring_CreateRes,
                  sstring_DeleteRes,
                  HammingIndep > >
       ( m_tests, m_gen,
         StatTest::Names({"Hamming Independence"}) ) );

  m_pvals.push_back(
    add< TestU01< swalk_Res,
                  swalk_CreateRes,
                  swalk_DeleteRes,
                  RandomWalk1 > >
       ( m_tests, m_gen,
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
  swrite_Basic = FALSE;         // Want screen no putput from TestU01

  //g_tid = 0;
  //bbattery_SmallCrush( m_gen.get() );

  using Pvals = StatTest::Pvals;
  using Psize = Pvals::size_type;
  using Tsize = TestContainer::size_type;

  Tsize ntest = m_tests.size();

  #ifdef _OPENMP
  #pragma omp parallel private(g_tid)
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
      Pvals pvals = m_tests[i]->run();
      // Evaluate test
      Psize npval = pvals.size();
      for (Psize p=0; p<npval; ++p) {
        if ((pvals[p] <= gofw_Suspectp) || (pvals[p] >= 1.0 - gofw_Suspectp)) {
          m_pvals[i][p] = pvals[p];
        }
      }
    }
  }

  // Output failed tests
  m_base.print.failed< StatTest >( "Failed tests", m_pvals, m_tests );

  m_base.print.endpart();
}
