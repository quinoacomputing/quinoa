//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.C
  \author    J. Bakosi
  \date      Thu 19 Jun 2014 11:00:35 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <SmallCrush.h>
#include <TestU01.h>

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::SmallCrush;

void
SmallCrush::addTests( std::vector< std::function< StatTest() > >& tests,
                      tk::ctr::RNGType rng,
                      CProxy_TestU01Suite& proxy )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
  // Select test stack
  const auto& stack = g_testStack.TestU01;

  // Find RNG properties
  auto gen = stack.rngprops(rng).ptr.get();

  static const long THOUSAND = 1000;
  static const long MILLION = THOUSAND * THOUSAND;
  static const long BILLION = THOUSAND * MILLION;

  // Marsaglia's BirthdaySpacings
  #ifdef USE_LONGLONG
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings"},
             1L, 5L * MILLION, 0, 1073741824L, 2, 1 );
  #else
  stack.add< TestU01< TestU01Props< tag::birthdayspacings, CProxy_TestU01Suite,
                                    sres_poisson, sres_createpoisson,
                                    sres_deletepoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"birthday spacings"},
             10l, million / 2, 0, 67108864l, 2, 1 );
  #endif

  // Knuth's Collision
  stack.add< TestU01< TestU01Props< tag::Collision, CProxy_TestU01Suite,
                                    sknuth_Res2, sknuth_CreateRes2,
                                    sknuth_DeleteRes2,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision"},
             1L, 5L * MILLION, 0, 65536L, 2 );

//   // Knuth's Gap
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( proxy, tests, id, gen, rng, {"Gap"},
//          Gap, 1L, MILLION / 5, 22, 0.0, 0.00390625 );
// 
//   // Knuth's Simple Poker
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( proxy, tests, id, gen, rng, {"Simplified Poker"},
//          SimpPoker, 1L, 2L * MILLION / 5, 24, 64, 64 );
// 
//   // Knuth's Coupon Collector
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( proxy, tests, id, gen, rng, {"Coupon Collector"},
//          CouponCollector, 1L, MILLION / 2, 26, 16 );
// 
//   // Knuth's Maximum-of-t
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( proxy, tests, id, gen, rng, {"Maximum-of-t", "Maximum-of-t Anderson-Darling"},
//          MaxOft, 1L, 2L * MILLION, 0, static_cast<int>(MILLION / 10), 6,
//                  gofw_Mean, gofw_Mean );
// 
//   // Weight Distribution
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( proxy, tests, id, gen, rng, {"Weight Distribution"},
//          WeightDistrib, 1L, MILLION / 5, 27, 256L, 0.0, 0.125 );
// 
//   // Marsaglia's Matrix Rank
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( proxy, tests, id, gen, rng, {"Matrix Rank"},
//          MatrixRank, 1L, 20L * THOUSAND, 20, 10, 60, 60 );
// 
//   // Hamming independence
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( proxy, tests, id, gen, rng, {"Hamming Independence"},
//          HammingIndep, 1L, MILLION/2, 20, 10, 300, 0 );

  // Random Walk 1
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
       ( proxy, tests, rng, gen, {"Random Walk 1 Stat H",
                                  "Random Walk 1 Stat M",
                                  "Random Walk 1 Stat J",
                                  "Random Walk 1 Stat R",
                                  "Random Walk 1 Stat C"},
         1L, static_cast<long>(MILLION), 0, 30, 150L, 150L );
}
