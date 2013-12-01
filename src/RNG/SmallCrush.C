//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Sat 30 Nov 2013 11:23:36 PM MST
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

extern std::vector< std::unique_ptr<tk::RNG> > g_rng;  //!< RNGs
extern Rsize g_rid;                                    //!< RNG id
extern std::vector< int > g_tid;                       //!< Thread ids

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

  // Instantiate all RNGs selected
  for (const auto& r : rngs) {
    if (r != ctr::RNGType::NO_RNG) {

      // Save name of random number generator
      std::string name( quinoa::ctr::RNG().name( r ) );

      // Instantiate random number generator
      g_rng.push_back( std::unique_ptr< tk::RNG >( m_base.rng[r]() ) );

      // Create associated TestU01 external generator
      char* const rngname = const_cast<char*>( name.c_str() );
      m_gen.push_back(
        Gen01Ptr( unif01_CreateExternGen01( rngname, MKLRNGUniform ) ) );

      // Add statistical tests to battery of RNG just instantiated
      m_tests.push_back( TestContainer() );
      m_pvals.push_back( std::vector< StatTest::Pvals >() );
      addTests( m_gen.back(), m_tests.back(), m_pvals.back() );
    }
  }
}

void
SmallCrush::addTests( const Gen01Ptr& gen,
                      TestContainer& tests,
                      std::vector< StatTest::Pvals >& pvals )
//******************************************************************************
// Add statistical tests to battery, count up total number of p-values
//! \author  J. Bakosi
//******************************************************************************
{
  pvals.push_back(
    add< TestU01< sres_Poisson,
                  sres_CreatePoisson,
                  sres_DeletePoisson,
                  BirthdaySpacings > >
       ( tests, gen, StatTest::Names({"Birthday Spacings"}) ) );
  

  pvals.push_back(
    add< TestU01< sknuth_Res2,
                  sknuth_CreateRes2,
                  sknuth_DeleteRes2,
                  Collision > >
       ( tests, gen, StatTest::Names({"Collision"}) ) );

  pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  Gap > >
       ( tests, gen, StatTest::Names({"Gap"}) ) );

  pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  SimpPoker > >
       ( tests, gen, StatTest::Names({"Simplified Poker"}) ) );

  pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  CouponCollector > >
       ( tests, gen, StatTest::Names({"Coupon Collector"}) ) );

  pvals.push_back(
    add< TestU01< sknuth_Res1,
                  sknuth_CreateRes1,
                  sknuth_DeleteRes1,
                  MaxOft > >
       ( tests, gen,
         StatTest::Names({"Maximum-of-t",
                          "Maximum-of-t Anderson-Darling"}) ) );

  pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  WeightDistrib > >
       ( tests, gen,
         StatTest::Names({"Weight Distribution"}) ) );

  pvals.push_back(
    add< TestU01< sres_Chi2,
                  sres_CreateChi2,
                  sres_DeleteChi2,
                  MatrixRank > >
       ( tests, gen, StatTest::Names({"Matrix Rank"}) ) );

  pvals.push_back(
    add< TestU01< sstring_Res,
                  sstring_CreateRes,
                  sstring_DeleteRes,
                  HammingIndep > >
       ( tests, gen,
         StatTest::Names({"Hamming Independence"}) ) );

  pvals.push_back(
    add< TestU01< swalk_Res,
                  swalk_CreateRes,
                  swalk_DeleteRes,
                  RandomWalk1 > >
       ( tests, gen,
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
//   const RNGTestPrint& print = m_base.print;
// 
//   print.part("Run SmallCrush");
//   swrite_Basic = FALSE;         // Want screen no putput from TestU01
// 
//   //g_tid = 0;
//   //bbattery_SmallCrush( m_gen.get() );
// 
//   using Pvals = StatTest::Pvals;
//   using Psize = Pvals::size_type;
//   using Tsize = TestContainer::size_type;
// 
//   Tsize ntest = m_tests.size();
// 
//   #ifdef _OPENMP
//   #pragma omp parallel private(g_tid)
//   #endif
//   {
//     #ifdef _OPENMP
//     g_tid = omp_get_thread_num();
//     #else
//     g_tid = 0;
//     #endif
// 
//     #ifdef _OPENMP
//     #pragma omp for
//     #endif
//     for (Tsize i=0; i<ntest; ++i) {
//       // Run test
//       Pvals pvals = m_tests[i]->run();
//       // Evaluate test
//       Psize npval = pvals.size();
//       for (Psize p=0; p<npval; ++p) {
//         if ((pvals[p] <= gofw_Suspectp) || (pvals[p] >= 1.0 - gofw_Suspectp)) {
//           m_pvals[i][p] = pvals[p];
//         }
//       }
//     }
//   }
// 
//   // Output failed tests
//   m_base.print.failed< StatTest >( "Failed tests", m_pvals, m_tests );
// 
//   m_base.print.endpart();
}

void
SmallCrush::print() const
//******************************************************************************
//  Print list of registered statistical tests
//! \author  J. Bakosi
//******************************************************************************
{
  // Output test names and number of p-values produced.
  // The first one suffices as all are the same.
  m_base.print.names< StatTest >( m_tests[0] );
}
