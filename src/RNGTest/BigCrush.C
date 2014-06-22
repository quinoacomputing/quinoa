//******************************************************************************
/*!
  \file      src/RNGTest/BigCrush.C
  \author    J. Bakosi
  \date      Mon 16 Jun 2014 11:21:41 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************

#include <BigCrush.h>
#include <TestU01.h>

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::BigCrush;

void
BigCrush::addTests( std::vector< std::function< StatTest() > >& tests,
                    tk::ctr::RNGType rng,
                    CProxy_TestU01Suite& proxy )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
//   // Marsaglia Serial Over, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Serial Over r=0"} ),
//        SerialOver, 1L, static_cast<long>(BILLION), 0, 256L, 3 );
// 
//   // Marsaglia Serial Over, r = 22
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Serial Over r=22"} ),
//        SerialOver, 1L, static_cast<long>(BILLION), 22, 256L, 3 );
// 
//   // Marsaglia Collision Over, t = 2, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=2 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 1024L * 1024L * 2L, 2 );
// 
//   // Marsaglia Collision Over, t = 2, r = 9
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=2 r=9"} ),
//        CollisionOver, 30L, 20L * MILLION, 9, 1024L * 1024L * 2L, 2 );
// 
//   // Marsaglia Collision Over, t = 3, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=3 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 1024L * 16L, 3 );
// 
//   // Marsaglia Collision Over, t = 3, r = 16
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=3 r=16"} ),
//        CollisionOver, 30L, 20L * MILLION, 16, 1024L * 16L, 3 );
// 
//   // Marsaglia Collision Over, t = 7, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=7 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 64L, 7 );
// 
//   // Marsaglia Collision Over, t = 7, r = 24
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=7 r=24"} ),
//        CollisionOver, 30L, 20L * MILLION, 24, 64L, 7 );
// 
//   // Marsaglia Collision Over, t = 14, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=14 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 8L, 14 );
// 
//   // Marsaglia Collision Over, t = 14, r = 27
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=14 r=27"} ),
//        CollisionOver, 30L, 20L * MILLION, 27, 8L, 14 );
// 
//   // Marsaglia Collision Over, t = 21, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=21 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 4L, 21 );
// 
//   // Marsaglia Collision Over, t = 21, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=21 r=0"} ),
//        CollisionOver, 30L, 20L * MILLION, 0, 4L, 21 );
// 
//   // Marsaglia Collision Over, t = 21, r = 28
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=21 r=28"} ),
//        CollisionOver, 30L, 20L * MILLION, 28, 4L, 21 );
// 
//   #ifdef USE_LONGLONG
// 
//   // Marsaglia's Birthday Spacings, t = 2, r = 0
//   #if LONG_MAX <= 2147483647L
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 250L, 4L * MILLION, 0, 1073741824L, 2, 1 );
//   #else // LONG_MAX
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 100L, 10L * MILLION, 0, 2147483648L, 2, 1 );
//   #endif // LONG_MAX
// 
//   // Marsaglia's Birthday Spacings, t = 3, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
//        BirthdaySpacings, 20L, 20L * MILLION, 0, 2097152L, 3, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 4, r = 14
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=14"} ),
//        BirthdaySpacings, 20L, 30L * MILLION, 14, 65536L, 4, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 7, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=0"} ),
//        BirthdaySpacings, 20L, 20L * MILLION, 0, 512L, 7, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 7, r = 7
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=7"} ),
//        BirthdaySpacings, 20L, 20L * MILLION, 7, 512L, 7, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 8, r = 14
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=14"} ),
//        BirthdaySpacings, 20L, 30L * MILLION, 14, 256L, 8, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 8, r = 22
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=22"} ),
//        BirthdaySpacings, 20L, 30L * MILLION, 22, 256L, 8, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 16, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=16 r=0"} ),
//        BirthdaySpacings, 20L, 30L * MILLION, 0, 16L, 16, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 16, r = 26
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=16 r=26"} ),
//        BirthdaySpacings, 20L, 30L * MILLION, 26, 16L, 16, 1 );
// 
//   #else // USE_LONGLONG
// 
//   // Marsaglia's Birthday Spacings, t = 2, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
//                          67108864L, 2, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 4, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
//                          1024L * 8L, 4, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 4, r = 16
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=16"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 16,
//                          1024L * 8L, 4, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=0"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
//                          16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 5
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=5"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 5,
//                          16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 10
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=10"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 10,
//                          16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 15
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=15"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 15,
//                          16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 20
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=20"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 20,
//                          16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 26
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=26"} ),
//        BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 26,
//                          16L, 13, 1 );
// 
//   #endif // USE_LONGLONG
// 
//   // Close Pairs, t = 3
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=3",
//                                    "Close Pairs mNP t=3",
//                                    "Close Pairs mNP1 t=3",
//                                    "Close Pairs mNP2 t=3",
//                                    "Close Pairs mNJumps t=3",
//                                    "Close Pairs mNP2S t=3"} ),
//        ClosePairs, 30L, 6L * MILLION, 0, 3, 0, 30, 1 );
// 
//   // Close Pairs, t = 5
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=5",
//                                    "Close Pairs mNP t=5",
//                                    "Close Pairs mNP1 t=5",
//                                    "Close Pairs mNP2 t=5",
//                                    "Close Pairs mNJumps t=5",
//                                    "Close Pairs mNP2S t=5"} ),
//        ClosePairs, 20L, 4L * MILLION, 0, 5, 0, 30, 1 );
// 
//   // Close Pairs, t = 9
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=9",
//                                    "Close Pairs mNP t=9",
//                                    "Close Pairs mNP1 t=9",
//                                    "Close Pairs mNP2 t=9",
//                                    "Close Pairs mNJumps t=9",
//                                    "Close Pairs mNP2S t=9"} ),
//        ClosePairs, 10L, 3L * MILLION, 0, 9, 0, 30, 1 );
// 
//   // Close Pairs, t = 16
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=16",
//                                    "Close Pairs mNP t=16",
//                                    "Close Pairs mNP1 t=16",
//                                    "Close Pairs mNP2 t=16",
//                                    "Close Pairs mNJumps t=16",
//                                    "Close Pairs mNP2S t=16"} ),
//        ClosePairs, 5L, 2L * MILLION, 0, 16, 0, 30, 1 );
// 
//   // Knuth's Simple Poker, d = 8, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=8 r=0"} ),
//          SimpPoker, 1L, 400L * MILLION, 0, 8, 8 );
// 
//   // Knuth's Simple Poker, d = 8, r = 27
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=8 r=27"} ),
//          SimpPoker, 1L, 400L * MILLION, 27, 8, 8 );
// 
//   // Knuth's Simple Poker, d = 32, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=32 r=0"} ),
//          SimpPoker, 1L, 100L * MILLION, 0, 32, 32 );
// 
//   // Knuth's Simple Poker, d = 32, r = 25
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=32 r=25"} ),
//          SimpPoker, 1L, 100L * MILLION, 25, 32, 32 );
// 
//   // Knuth's Coupon Collector, d = 8, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=8 r=0"} ),
//          CouponCollector, 1L, 200L * MILLION, 0, 8 );
// 
//   // Knuth's Coupon Collector, d = 8, r = 10
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=8 r=10"} ),
//          CouponCollector, 1L, 200L * MILLION, 10, 8 );
// 
//   // Knuth's Coupon Collector, d = 8, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=8 r=20"} ),
//          CouponCollector, 1L, 200L * MILLION, 20, 8 );
// 
//   // Knuth's Coupon Collector, d = 8, r = 27
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=8 r=27"} ),
//          CouponCollector, 1L, 200L * MILLION, 27, 8 );
// 
//   // Knuth's Gap, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=0"} ),
//          Gap, 1L, static_cast<long>(BILLION/2), 0, 0.0, 1.0/16.0 );
// 
//   // Knuth's Gap, r = 25
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=25"} ),
//          Gap, 1L, 300L * MILLION, 25, 0.0, 1.0/32.0 );
// 
//   // Knuth's Gap, r = 0, n = .5e+9
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=0 n=.5B"} ),
//          Gap, 1L, static_cast<long>(BILLION / 10), 0, 0.0, 1.0/128.0 );
// 
//   // Knuth's Gap, r = 20, n = 10e+6
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=20 n=10M"} ),
//          Gap, 1L, 10L * MILLION, 20, 0.0, 1.0/1024.0 );
// 
//   // Knuth's Run, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Run r=0"} ),
//          Run, 5L, static_cast<long>(BILLION), 0, 0 );
// 
//   // Knuth's Run, r = 15
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Run r=15"} ),
//          Run, 10L, static_cast<long>(BILLION), 15, 1 );
// 
//   // Knuth's Permutation, t = 3
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation t=3"} ),
//          Permutation, 1L, static_cast<long>(BILLION), 5, 3 );
// 
//   // Knuth's Permutation, t = 5
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation t=5"} ),
//          Permutation, 1L, static_cast<long>(BILLION), 5, 5 );
// 
//   // Knuth's Permutation, t = 7
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation t=7"} ),
//          Permutation, 1L, static_cast<long>(BILLION / 2), 5, 7 );
// 
//   // Knuth's Permutation, t = 10
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation t=10"} ),
//          Permutation, 1L, static_cast<long>(BILLION / 2), 10, 10 );
// 
//   // Knuth's Collision with permutations, r = 0
//   add< TestU01< sknuth_Res2, sknuth_CreateRes2, sknuth_DeleteRes2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Collision w. Permutations r=0"} ),
//          CollisionPermut, 20L, 20L * MILLION, 0, 14 );
// 
//   // Knuth's Collision with permutations, r = 10
//   add< TestU01< sknuth_Res2, sknuth_CreateRes2, sknuth_DeleteRes2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Collision w. Permutations r=10"} ),
//          CollisionPermut, 20L, 20L * MILLION, 10, 14 );
// 
//   // Knuth's Maximum-of-t, t = 8
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=8",
//                                      "Maximum-of-t Anderson-Darling t=8"} ),
//          MaxOft, 40L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 8,
//                  gofw_Sum, gofw_AD );
// 
//   // Knuth's Maximum-of-t, t = 16
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=16",
//                                          "Maximum-of-t Anderson-Darling t=16"} ),
//          MaxOft, 30L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 16,
//                  gofw_Sum, gofw_AD );
// 
//   // Knuth's Maximum-of-t, t = 24
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=24",
//                                          "Maximum-of-t Anderson-Darling t=24"} ),
//          MaxOft, 20L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 24,
//                  gofw_Mean, gofw_Mean );
// 
//   // Knuth's Maximum-of-t, t = 32
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=32",
//                                          "Maximum-of-t Anderson-Darling t=32"} ),
//          MaxOft, 20L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 32,
//                  gofw_Mean, gofw_Mean );
// 
//   // Sample Products, t = 8
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Products t=8"} ),
//          SampleProd, 4L, 10L * MILLION, 0, 8 );
// 
//   // Sample Products, t = 16
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Products t=16"} ),
//          SampleProd, 20L, 10L * MILLION, 0, 16 );
// 
//   // Sample Products, t = 24
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Products t=24"} ),
//          SampleProd, 20L, 10L * MILLION, 0, 24 );
// 
//   // Sample Mean, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Mean r=0"} ),
//          SampleMean, 20L * MILLION, 30L, 0 );
// 
//   // Sample Mean, r = 10
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Mean r=10"} ),
//          SampleMean, 20L * MILLION, 30L, 10 );
// 
//   // Sample Autorrelation, k = 1
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Autorrelation k=1"} ),
//          SampleCorr, 1L, 2L * BILLION, 0, 1 );
// 
//   // Sample Autorrelation, k = 2
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Autorrelation k=2"} ),
//          SampleCorr, 1L, 2L * BILLION, 0, 2 );
// 
//   // Maurer's "universal" test, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Appearance Spacings r=0"} ),
//          AppearanceSpacings, 1L, 10L * MILLION, 400L * MILLION, 0, 30, 15 );
// 
//   // Maurer's "universal" test, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Appearance Spacings r=0"} ),
//          AppearanceSpacings, 1L, 10L * MILLION, static_cast<long>(BILLION), 0,
//                              3, 15 );
// 
//   // Maurer's "universal" test, r = 27
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Appearance Spacings r=27"} ),
//          AppearanceSpacings, 1L, 10L * MILLION, static_cast<long>(BILLION), 27,
//                              3, 15 );
// 
//   // Weight Distribution, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=0"} ),
//          WeightDistrib, 1L, 20L * MILLION, 0, 256L, 0.0, 0.25 );
// 
//   // Weight Distribution, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=20"} ),
//          WeightDistrib, 1L, 20L * MILLION, 20, 256L, 0.0, 0.25 );
// 
//   // Weight Distribution, r = 28
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=28"} ),
//          WeightDistrib, 1L, 20L * MILLION, 28, 256L, 0.0, 0.25 );
// 
//   // Weight Distribution, r = 0, beta = 1/16
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=0 beta=1/16"} ),
//          WeightDistrib, 1L, 20L * MILLION, 0, 256L, 0.0, 0.0625 );
// 
//   // Weight Distribution, r = 10
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=10"} ),
//          WeightDistrib, 1L, 20L * MILLION, 10, 256L, 0.0, 0.0625 );
// 
//   // Weight Distribution, r = 26
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=26"} ),
//          WeightDistrib, 1L, 20L * MILLION, 26, 256L, 0.0, 0.0625 );
// 
//   // Sum Collector
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sum Collector"} ),
//          SumCollector, 1L, 500L * MILLION, 0, 10.0 );
// 
//   // Marsaglia's Matrix Rank, L = 30, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=30 r=0"} ),
//          MatrixRank, 10L, static_cast<long>(MILLION), 0, 5, 30, 30 );
// 
//   // Marsaglia's Matrix Rank, L = 30, r = 25
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=30 r=25"} ),
//          MatrixRank, 10L, static_cast<long>(MILLION), 25, 5, 30, 30 );
// 
//   // Marsaglia's Matrix Rank, L = 1000, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=1000 r=0"} ),
//          MatrixRank, 1L, 5L * THOUSAND, 0, 4, 1000, 1000 );
// 
//   // Marsaglia's Matrix Rank, L = 1000, r = 26
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=1000 r=26"} ),
//          MatrixRank, 1L, 5L * THOUSAND, 26, 4, 1000, 1000 );
// 
//   // Marsaglia's Matrix Rank, L = 5000, r = 15
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=5000 r=15"} ),
//        MatrixRank, 1L, 80L, 15, 15, 5000, 5000 );
// 
//   // Marsaglia's Matrix Rank, L = 5000, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank L=5000 r=0"} ),
//        MatrixRank, 1L, 80L, 0, 30, 5000, 5000 );
// 
//   // Marsaglia's Modified Savir
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Savir2"} ),
//        Savir2, 10L, 10L * MILLION, 10, 1024L * 1024L, 30 );
// 
//   // Marsaglia's greatest common divisor
//   add< TestU01< smarsa_Res2, smarsa_CreateRes2, smarsa_DeleteRes2,
//               long, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"GCD"} ),
//        GCD, 10L, 50L * MILLION, 0, 30 );
// 
//   // Random Walk 1, L = 50, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=50 r=0",
//                                          "Random Walk 1 Stat M L=50 r=0",
//                                          "Random Walk 1 Stat J L=50 r=0",
//                                          "Random Walk 1 Stat R L=50 r=0",
//                                          "Random Walk 1 Stat C L=50 r=0"} ),
//          RandomWalk1, 1L, 100L * MILLION, 0, 5, 50L, 50L );
// 
//   // Random Walk 1, L = 50, r = 25
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=50 r=25",
//                                          "Random Walk 1 Stat M L=50 r=25",
//                                          "Random Walk 1 Stat J L=50 r=25",
//                                          "Random Walk 1 Stat R L=50 r=25",
//                                          "Random Walk 1 Stat C L=50 r=25"} ),
//          RandomWalk1, 1L, 100L * MILLION, 25, 5, 50L, 50L );
// 
//   // Random Walk 1, L = 1000, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=1000 r=0",
//                                          "Random Walk 1 Stat M L=1000 r=0",
//                                          "Random Walk 1 Stat J L=1000 r=0",
//                                          "Random Walk 1 Stat R L=1000 r=0",
//                                          "Random Walk 1 Stat C L=1000 r=0"} ),
//          RandomWalk1, 1L, 10L * MILLION, 0, 10, 1000L, 1000L );
// 
//   // Random Walk 1, L = 1000, r = 20
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=1000 r=20",
//                                          "Random Walk 1 Stat M L=1000 r=20",
//                                          "Random Walk 1 Stat J L=1000 r=20",
//                                          "Random Walk 1 Stat R L=1000 r=20",
//                                          "Random Walk 1 Stat C L=1000 r=20"} ),
//          RandomWalk1, 1L, 10L * MILLION, 20, 10, 1000L, 1000L );
// 
//   // Random Walk 1, L = 10000, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=10000 r=0",
//                                          "Random Walk 1 Stat M L=10000 r=0",
//                                          "Random Walk 1 Stat J L=10000 r=0",
//                                          "Random Walk 1 Stat R L=10000 r=0",
//                                          "Random Walk 1 Stat C L=10000 r=0"} ),
//          RandomWalk1, 1L, static_cast<long>(MILLION), 0, 15, 10000L, 10000L );
// 
//   // Random Walk 1, L = 10000, r = 15
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=10000 r=15",
//                                          "Random Walk 1 Stat M L=10000 r=15",
//                                          "Random Walk 1 Stat J L=10000 r=15",
//                                          "Random Walk 1 Stat R L=10000 r=15",
//                                          "Random Walk 1 Stat C L=10000 r=15"} ),
//          RandomWalk1, 1L, static_cast<long>(MILLION), 15, 15, 10000L, 10000L );
// 
//   // Linear Complexity, r = 0
//   add< TestU01< scomp_Res, scomp_CreateRes, scomp_DeleteRes,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Linear Complexity Jump r=0",
//                                          "Linear Complexity Size r=0"} ),
//          LinearComp, 1L, 400L * THOUSAND + 20, 0, 1 );
// 
//   // Linear Complexity, r = 29
//   add< TestU01< scomp_Res, scomp_CreateRes, scomp_DeleteRes,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Linear Complexity Jump r=29",
//                                          "Linear Complexity Size r=29"} ),
//          LinearComp, 1L, 400L * THOUSAND + 20, 29, 1 );
// 
//   // Lempel-Ziv Compressibility, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Lempel-Ziv Compressibility r=0"} ),
//          LempelZiv, 10L, 27, 0, 30 );
// 
//   // Lempel-Ziv Compressibility, r = 15
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Lempel-Ziv Compressibility r=15"} ),
//          LempelZiv, 10L, 27, 15, 15 );
// 
//   // Fourier3, r = 0
//   add< TestU01< sspectral_Res, sspectral_CreateRes, sspectral_DeleteRes,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Fourier 3 r=0"} ),
//          Fourier3, 100L * THOUSAND, 14, 0, 3 );
// 
//   // Fourier3, r = 20
//   add< TestU01< sspectral_Res, sspectral_CreateRes, sspectral_DeleteRes,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Fourier 3 r=27"} ),
//          Fourier3, 100L * THOUSAND, 14, 27, 3 );
// 
//   // Longes Heat Run, r = 0
//   add< TestU01< sstring_Res2, sstring_CreateRes2, sstring_DeleteRes2,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Longest Head Run Chi r=0",
//                                          "Longest Head Run Disc r=0"} ),
//          LongestHeadRun, 1L, 1000L, 0, 3, 20L + 10L * MILLION );
// 
//   // Longes Heat Run, r = 27
//   add< TestU01< sstring_Res2, sstring_CreateRes2, sstring_DeleteRes2,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Longest Head Run Chi r=27",
//                                          "Longest Head Run Disc r=27"} ),
//          LongestHeadRun, 1L, 1000L, 27, 3, 20L + 10L * MILLION );
// 
//   // Periods In Strings, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Periods In Strings r=0"} ),
//          PeriodsInStrings, 10L, BILLION / 2, 0, 10 );
// 
//   // Periods In Strings, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Periods In Strings r=20"} ),
//          PeriodsInStrings, 10L, BILLION / 2, 20, 10 );
// 
//   // Hamming Weight 2, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Weight 2 r=0"} ),
//          HammingWeight2, 10L, static_cast<long>(BILLION), 0, 3,
//                          static_cast<long>(MILLION) );
// 
//   // Hamming Weight 2, r = 27
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Weight 2 r=27"} ),
//          HammingWeight2, 10L, static_cast<long>(BILLION), 27, 3,
//                          static_cast<long>(MILLION) );
// 
//   // Hamming Correlation, L = 30
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=30"} ),
//          HammingCorr, 1L, static_cast<long>(BILLION), 10, 10, 30 );
// 
//   // Hamming Correlation, L = 300
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=300"} ),
//          HammingCorr, 1L, 100L * MILLION, 10, 10, 10 * 30 );
// 
//   // Hamming Correlation, L = 1200
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=1200"} ),
//          HammingCorr, 1L, 100L * MILLION, 10, 10, 40 * 30 );
// 
//   // Hamming independence, L = 30, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=30 r=0"} ),
//          HammingIndep, 10L, 30L * MILLION, 0, 3, 30, 0 );
// 
//   // Hamming independence, L = 30, r = 27
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=30 r=27"} ),
//          HammingIndep, 10L, 30L * MILLION, 27, 3, 30, 0 );
// 
//   // Hamming independence, L = 300, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=300 r=0"} ),
//          HammingIndep, 1L, 30L * MILLION, 0, 4, 10 * 30, 0 );
// 
//   // Hamming independence, L = 300, r = 26
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=300 r=26"} ),
//          HammingIndep, 1L, 30L * MILLION, 26, 4, 10 * 30, 0 );
// 
//   // Hamming independence, L = 1200, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=1200 r=0"} ),
//          HammingIndep, 1L, 10L * MILLION, 0, 5, 40 * 30, 0 );
// 
//   // Hamming independence, L = 1200, r = 25
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=1200 r=25"} ),
//          HammingIndep, 1L, 10L * MILLION, 25, 5, 40 * 30, 0 );
// 
//   // String Run, r = 0
//   add< TestU01< sstring_Res3, sstring_CreateRes3, sstring_DeleteRes3,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"String Run NRuns r=0",
//                                          "String Run NBits r=0"} ),
//          StringRun, 1L, 2L * BILLION, 0, 3 );
// 
//   // String Run, r = 27
//   add< TestU01< sstring_Res3, sstring_CreateRes3, sstring_DeleteRes3,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"String Run NRuns r=27",
//                                          "String Run NBits r=27"} ),
//          StringRun, 1L, 2L * BILLION, 27, 3 );
// 
//   // Autocorrelation, d = 1, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=1 r=0"} ),
//          AutoCor, 10L, 30L + BILLION, 0, 3, 1 );
// 
//   // Autocorrelation, d = 3, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=3 r=0"} ),
//          AutoCor, 10L, 30L + BILLION, 0, 3, 3 );
// 
//   // Autocorrelation, d = 1, r =27
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=1 r=27"} ),
//          AutoCor, 10L, 30L + BILLION, 27, 3, 1 );
// 
//   // Autocorrelation, d = 3, r = 27
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=3 r=27"} ),
//          AutoCor, 10L, 30L + BILLION, 27, 3, 3 );
}
