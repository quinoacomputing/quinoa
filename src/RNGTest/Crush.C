//******************************************************************************
/*!
  \file      src/RNGTest/Crush.C
  \author    J. Bakosi
  \date      Sat 07 Jun 2014 06:36:05 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Crush battery
  \details   Crush battery
*/
//******************************************************************************

#include <Crush.h>

using rngtest::Crush;

void
Crush::addTests( std::vector< StatTest >& tests, tk::ctr::RNGType r )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
//   // Marsaglia's Serial Over, t = 2
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Serial Over t=2"} ),
//        SerialOver, 1L, 500L * MILLION, 0, 4096L, 2 );
// 
//   // Marsaglia's Serial Over, t = 4
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Serial Over t=4"} ),
//        SerialOver, 1L, 300L * MILLION, 0, 64L, 4 );
// 
//   // Marsaglia's Collision Over, t = 2, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=2 r=0"} ),
//        CollisionOver, 10L, 10L * MILLION, 0, 1024L * 1024, 2 );
// 
//   // Marsaglia's Collision Over, t = 2, r = 10
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=2 r=10"} ),
//        CollisionOver, 10L, 10L * MILLION, 10, 1024L * 1024, 2 );
// 
//   // Marsaglia's Collision Over, t = 4, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=4 r=0"} ),
//        CollisionOver, 10L, 10L * MILLION, 0, 1024L, 4 );
// 
//   // Marsaglia's Collision Over, t = 4, r = 20
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=4 r=20"} ),
//        CollisionOver, 10L, 10L * MILLION, 20, 1024L, 4 );
// 
//   // Marsaglia's Collision Over, t = 8, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=8 r=0"} ),
//        CollisionOver, 10L, 10L * MILLION, 0, 32L, 8 );
// 
//   // Marsaglia's Collision Over, t = 8, r = 25
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=8 r=25"} ),
//        CollisionOver, 10L, 10L * MILLION, 25, 32L, 8 );
// 
//   // Marsaglia's Collision Over, t = 20, r = 0
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=20 r=0"} ),
//        CollisionOver, 10L, 10L * MILLION, 0, 4L, 20 );
// 
//   // Marsaglia's Collision Over, t = 20, r = 28
//   add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
//                 long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Collision Over t=20 r=28"} ),
//        CollisionOver, 10L, 10L * MILLION, 28, 4L, 20 );
// 
//   #ifdef USE_LONGLONG
// 
//   // Marsaglia's Birthday Spacings, t = 2, r = 0
//   #if LONG_MAX <= 2147483647L
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 10L, 10L * MILLION, 0, 1073741824L, 2, 1 );
//   #else // LONG_MAX
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 5L, 20L * MILLION, 0, 2L*1073741824L, 2, 1 );
//   #endif // LONG_MAX
// 
//   // Marsaglia's Birthday Spacings, t = 3, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
//        BirthdaySpacings, 5L, 20L * MILLION, 0, 2097152L, 3, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 4, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
//        BirthdaySpacings, 5L, 20L * MILLION, 0, 65536L, 4, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 7, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=0"} ),
//        BirthdaySpacings, 3L, 20L * MILLION, 0, 512L, 7, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 7, r = 7
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=7"} ),
//        BirthdaySpacings, 3L, 20L * MILLION, 7, 512L, 7, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 8, r = 14
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=14"} ),
//        BirthdaySpacings, 3L, 20L * MILLION, 14, 256L, 8, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 8, r = 22
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=22"} ),
//        BirthdaySpacings, 3L, 20L * MILLION, 22, 256L, 8, 1 );
// 
//   #else // USE_LONGLONG
// 
//   // Marsaglia's Birthday Spacings, t = 2, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
//        BirthdaySpacings, 200L, 4L * MILLION / 10, 0, 67108864L, 2, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 3, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
//        BirthdaySpacings, 100L, 4L * MILLION / 10, 0, 131072L, 3, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 4, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
//        BirthdaySpacings, 200L, 4L * MILLION / 10, 0, 1024L * 8, 4, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 0
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=0"} ),
//        BirthdaySpacings, 100L, 4L * MILLION / 10, 0, 16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 10
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=10"} ),
//        BirthdaySpacings, 100L, 4L * MILLION / 10, 10, 16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 20
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=20"} ),
//        BirthdaySpacings, 100L, 4L * MILLION / 10, 20, 16L, 13, 1 );
// 
//   // Marsaglia's Birthday Spacings, t = 13, r = 26
//   add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
//                 long, long, int, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=26"} ),
//        BirthdaySpacings, 100L, 4L * MILLION / 10, 26, 16L, 13, 1 );
// 
//   #endif // USE_LONGLONG
// 
//   // Close Pairs, t = 2
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=2",
//                                        "Close Pairs mNP t=2",
//                                        "Close Pairs mNP1 t=2",
//                                        "Close Pairs mNP2 t=2",
//                                        "Close Pairs mNJumps t=2"} ),
//        ClosePairs, 10L, 2L * MILLION, 0, 2, 0, 30, 0 );
// 
//   // Close Pairs, t = 3
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=3",
//                                        "Close Pairs mNP t=3",
//                                        "Close Pairs mNP1 t=3",
//                                        "Close Pairs mNP2 t=3",
//                                        "Close Pairs mNJumps t=3",
//                                        "Close Pairs mNP2S t=3"} ),
//        ClosePairs, 10L, 2L * MILLION, 0, 3, 0, 30, 1 );
// 
//   // Close Pairs, t = 7
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs NP t=3",
//                                        "Close Pairs mNP t=3",
//                                        "Close Pairs mNP1 t=3",
//                                        "Close Pairs mNP2 t=3",
//                                        "Close Pairs mNJumps t=3",
//                                        "Close Pairs mNP2S t=3"} ),
//        ClosePairs, 5L, 2L * MILLION, 0, 7, 0, 30, 1 );
// 
//   // Close Pairs Bit Match, t = 2
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs Bit Match t=2"} ),
//        ClosePairsBitMatch, 4L, 4L * MILLION, 0, 2 );
// 
//   // Close Pairs Bit Match, t = 4
//   add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
//                 long, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Close Pairs Bit Match t=4"} ),
//        ClosePairsBitMatch, 2L, 4L * MILLION, 0, 4 );
// 
//   // Knuth's Simple Poker, d = 16, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=16 r=0"} ),
//          SimpPoker, 1L, 40L * MILLION, 0, 16, 16 );
// 
//   // Knuth's Simple Poker, d = 16, r = 26
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=16 r=26"} ),
//          SimpPoker, 1L, 40L * MILLION, 26, 16, 16 );
// 
//   // Knuth's Simple Poker, d = 64, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=64 r=0"} ),
//          SimpPoker, 1L, 10L * MILLION, 0, 64, 64 );
// 
//   // Knuth's Simple Poker, d = 64, r = 24
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Simplified Poker d=64 r=24"} ),
//          SimpPoker, 1L, 10L * MILLION, 24, 64, 64 );
// 
//   // Knuth's Coupon Collector, d = 4, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=4 r=0"} ),
//          CouponCollector, 1L, 40L * MILLION, 0, 4 );
// 
//   // Knuth's Coupon Collector, d = 4, r = 28
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=4 r=28"} ),
//          CouponCollector, 1L, 40L * MILLION, 28, 4 );
// 
//   // Knuth's Coupon Collector, d = 16, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=16 r=0"} ),
//          CouponCollector, 1L, 10L * MILLION, 0, 16 );
// 
//   // Knuth's Coupon Collector, d = 16, r = 26
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Coupon Collector d=16 r=26"} ),
//          CouponCollector, 1L, 10L * MILLION, 26, 16 );
// 
//   // Knuth's Gap, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=0"} ),
//          Gap, 1L, 100L * MILLION, 0, 0.0, 0.125 );
// 
//   // Knuth's Gap, r = 27
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=27"} ),
//          Gap, 1L, 100L * MILLION, 27, 0.0, 0.125 );
// 
//   // Knuth's Gap, r = 0, n = 5e+6
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=0 n=5M"} ),
//          Gap, 1L, 5L * MILLION, 0, 0.0, 1.0/256.0 );
// 
//   // Knuth's Gap, r = 22, n = 5e+6
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Gap r=22 n=5M"} ),
//          Gap, 1L, 5L * MILLION, 22, 0.0, 1.0/256.0 );
// 
//   // Knuth's Run, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Run r=0"} ),
//          Run, 1L, 500L * MILLION, 0, 1 );
// 
//   // Knuth's Run, r = 15
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Run r=15"} ),
//          Run, 1L, 500L * MILLION, 15, 0 );
// 
//   // Knuth's Permutation, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation r=0"} ),
//          Permutation, 1L, 50L * MILLION, 0, 10 );
// 
//   // Knuth's Permutation, r = 15
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Permutation r=15"} ),
//          Permutation, 1L, 50L * MILLION, 15, 10 );
// 
//   // Knuth's Collision with permutations, r = 0
//   add< TestU01< sknuth_Res2, sknuth_CreateRes2, sknuth_DeleteRes2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Collision w. Permutations r=0"} ),
//          CollisionPermut, 5L, 10L * MILLION, 0, 13 );
// 
//   // Knuth's Collision with permutations, r = 15
//   add< TestU01< sknuth_Res2, sknuth_CreateRes2, sknuth_DeleteRes2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Collision w. Permutations r=15"} ),
//          CollisionPermut, 5L, 10L * MILLION, 15, 13 );
// 
//   // Knuth's Maximum-of-t, t = 5
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=5",
//                                          "Maximum-of-t Anderson-Darling t=5"} ),
//          MaxOft, 10L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 5,
//                  gofw_Sum, gofw_AD );
// 
//   // Knuth's Maximum-of-t, t = 10
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=10",
//                                          "Maximum-of-t Anderson-Darling t=10"} ),
//          MaxOft, 5L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 10,
//                  gofw_Sum, gofw_AD );
// 
//   // Knuth's Maximum-of-t, t = 20
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=20",
//                                          "Maximum-of-t Anderson-Darling t=20"} ),
//          MaxOft, 1L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 20,
//                  gofw_Mean, gofw_Mean );
// 
//   // Knuth's Maximum-of-t, t = 30
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( tests, id, gen, rng, StatTest::Names( {"Maximum-of-t t=30",
//                                          "Maximum-of-t Anderson-Darling t=30"} ),
//          MaxOft, 1L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 30,
//                  gofw_Mean, gofw_Mean );
// 
//   // Sample Products, t = 10
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Products t=10"} ),
//          SampleProd, 1L, 10L * MILLION, 0, 10 );
// 
//   // Sample Products, t = 30
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Products t=30"} ),
//          SampleProd, 1L, 10L * MILLION, 0, 30 );
// 
//   // Sample Mean
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Mean"} ),
//          SampleMean, 10L * MILLION, 20L, 0 );
// 
//   // Sample Autorrelation
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sample Autorrelation"} ),
//          SampleCorr, 1L, 500L * MILLION, 0, 1 );
// 
//   // Maurer's "universal" test, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Appearance Spacings r=0"} ),
//          AppearanceSpacings, 1L, 10L * MILLION, 400L * MILLION, 0, 30, 15 );
// 
//   // Maurer's "universal" test, r = 20
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Appearance Spacings r=20"} ),
//          AppearanceSpacings, 1L, 10L * MILLION, 100L * MILLION, 20, 10, 15 );
// 
//   // Weight Distribution, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=0"} ),
//          WeightDistrib, 1L, 2L * MILLION, 0, 256L, 0.0, 0.125 );
// 
//   // Weight Distribution, r = 8
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=8"} ),
//          WeightDistrib, 1L, 2L * MILLION, 8, 256L, 0.0, 0.125 );
// 
//   // Weight Distribution, r = 16
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=16"} ),
//          WeightDistrib, 1L, 2L * MILLION, 16, 256L, 0.0, 0.125 );
// 
//   // Weight Distribution, r = 24
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Weight Distribution r=24"} ),
//          WeightDistrib, 1L, 2L * MILLION, 24, 256L, 0.0, 0.125 );
// 
//   // Sum Collector
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double > >
//        ( tests, id, gen, rng, StatTest::Names( {"Sum Collector"} ),
//          SumCollector, 1L, 20L * MILLION, 0, 10.0 );
// 
//   // Marsaglia's Matrix Rank, 60 x 60, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 60x60 r=0"} ),
//          MatrixRank, 1L, static_cast<long>(MILLION), 0, 30, 2 * 30, 2 * 30 );
// 
//   // Marsaglia's Matrix Rank, 60 x 60, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 60x60 r=20"} ),
//          MatrixRank, 1L, static_cast<long>(MILLION), 20, 10, 2 * 30, 2 * 30 );
// 
//   // Marsaglia's Matrix Rank, 300 x 300, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 300x300 r=0"} ),
//          MatrixRank, 1L, 50L * THOUSAND, 0, 30, 10 * 30, 10 * 30 );
// 
//   // Marsaglia's Matrix Rank, 300 x 300, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 300x300 r=20"} ),
//          MatrixRank, 1L, 50L * THOUSAND, 20, 10, 10 * 30, 10 * 30 );
// 
//   // Marsaglia's Matrix Rank, 1200 x 1200, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 1200x1200 r=0"} ),
//        MatrixRank, 1L, 2L * THOUSAND, 0, 30, 40 * 30, 40 * 30 );
// 
//   // Marsaglia's Matrix Rank, 1200 x 1200, r = 20
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, int, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Matrix Rank 1200x1200 r=20"} ),
//        MatrixRank, 1L, 2L * THOUSAND, 20, 10, 40 * 30, 40 * 30 );
// 
//   // Marsaglia's Modified Savir
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//               long, long, int, long, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"Savir2"} ),
//        Savir2, 1L, 20L * MILLION, 0, 1024L * 1024L, 30 );
// 
//   // Marsaglia's greatest common divisor, r = 0
//   add< TestU01< smarsa_Res2, smarsa_CreateRes2, smarsa_DeleteRes2,
//               long, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"GCD r=0"} ),
//        GCD, 1L, 100L * MILLION, 0, 30 );
// 
//   // Marsaglia's greatest common divisor, r = 10
//   add< TestU01< smarsa_Res2, smarsa_CreateRes2, smarsa_DeleteRes2,
//               long, long, int, int > >
//      ( tests, id, gen, rng, StatTest::Names( {"GCD r=10"} ),
//        GCD, 1L, 40L * MILLION, 10, 20 );
// 
//   // Random Walk 1, L = 90, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=90 r=0",
//                                          "Random Walk 1 Stat M L=90 r=0",
//                                          "Random Walk 1 Stat J L=90 r=0",
//                                          "Random Walk 1 Stat R L=90 r=0",
//                                          "Random Walk 1 Stat C L=90 r=0"} ),
//          RandomWalk1, 1L, 50L * MILLION, 0, 30, 90L, 90L );
// 
//   // Random Walk 1, L = 90, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=90 r=20",
//                                          "Random Walk 1 Stat M L=90 r=20",
//                                          "Random Walk 1 Stat J L=90 r=20",
//                                          "Random Walk 1 Stat R L=90 r=20",
//                                          "Random Walk 1 Stat C L=90 r=20"} ),
//          RandomWalk1, 1L, 10L * MILLION, 20, 10, 90L, 90L );
// 
//   // Random Walk 1, L = 1000, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=1000 r=0",
//                                          "Random Walk 1 Stat M L=1000 r=0",
//                                          "Random Walk 1 Stat J L=1000 r=0",
//                                          "Random Walk 1 Stat R L=1000 r=0",
//                                          "Random Walk 1 Stat C L=1000 r=0"} ),
//          RandomWalk1, 1L, 5L * MILLION, 0, 30, 1000L, 1000L );
// 
//   // Random Walk 1, L = 1000, r = 20
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=1000 r=20",
//                                          "Random Walk 1 Stat M L=1000 r=20",
//                                          "Random Walk 1 Stat J L=1000 r=20",
//                                          "Random Walk 1 Stat R L=1000 r=20",
//                                          "Random Walk 1 Stat C L=1000 r=20"} ),
//          RandomWalk1, 1L, static_cast<long>(MILLION), 20, 10, 1000L, 1000L );
// 
//   // Random Walk 1, L = 10000, r = 0
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=10000 r=0",
//                                          "Random Walk 1 Stat M L=10000 r=0",
//                                          "Random Walk 1 Stat J L=10000 r=0",
//                                          "Random Walk 1 Stat R L=10000 r=0",
//                                          "Random Walk 1 Stat C L=10000 r=0"} ),
//          RandomWalk1, 1L, MILLION / 2, 0, 30, 10000L, 10000L );
// 
//   // Random Walk 1, L = 10000, r = 20
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H L=10000 r=20",
//                                          "Random Walk 1 Stat M L=10000 r=20",
//                                          "Random Walk 1 Stat J L=10000 r=20",
//                                          "Random Walk 1 Stat R L=10000 r=20",
//                                          "Random Walk 1 Stat C L=10000 r=20"} ),
//          RandomWalk1, 1L, MILLION / 10, 20, 10, 10000L, 10000L );
// 
//   // Linear Complexity, r = 0
//   add< TestU01< scomp_Res, scomp_CreateRes, scomp_DeleteRes,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Linear Complexity Jump r=0",
//                                      "Linear Complexity Size r=0"} ),
//          LinearComp, 1L, 120L * THOUSAND, 0, 1 );
// 
//   // Linear Complexity, r = 29
//   add< TestU01< scomp_Res, scomp_CreateRes, scomp_DeleteRes,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Linear Complexity Jump r=29",
//                                          "Linear Complexity Size r=29"} ),
//          LinearComp, 1L, 120L * THOUSAND, 29, 1 );
// 
//   // Lempel-Ziv Compressibility
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Lempel-Ziv Compressibility"} ),
//          LempelZiv, 10L, 25, 0, 30 );
// 
//   // Fourier3, r = 0
//   add< TestU01< sspectral_Res, sspectral_CreateRes, sspectral_DeleteRes,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Fourier 3 r=0"} ),
//          Fourier3, 50L * THOUSAND, 14, 0, 30 );
// 
//   // Fourier3, r = 20
//   add< TestU01< sspectral_Res, sspectral_CreateRes, sspectral_DeleteRes,
//                 long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Fourier 3 r=20"} ),
//          Fourier3, 50L * THOUSAND, 14, 20, 10 );
// 
//   // Longes Heat Run, r = 0
//   add< TestU01< sstring_Res2, sstring_CreateRes2, sstring_DeleteRes2,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Longest Head Run Chi r=0",
//                                          "Longest Head Run Disc r=0"} ),
//          LongestHeadRun, 1L, 1000L, 0, 30, 20L + 10L * MILLION );
// 
//   // Longes Heat Run, r = 20
//   add< TestU01< sstring_Res2, sstring_CreateRes2, sstring_DeleteRes2,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Longest Head Run Chi r=20",
//                                          "Longest Head Run Disc r=20"} ),
//          LongestHeadRun, 1L, 300L, 20, 10, 20L + 10L * MILLION );
// 
//   // Periods In Strings, r = 0
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Periods In Strings r=0"} ),
//          PeriodsInStrings, 1L, 300L * MILLION, 0, 30 );
// 
//   // Periods In Strings, r = 15
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Periods In Strings r=15"} ),
//          PeriodsInStrings, 1L, 300L * MILLION, 15, 15 );
// 
//   // Hamming Weight 2, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Weight 2 r=0"} ),
//          HammingWeight2, 100L, 100L * MILLION, 0, 30,
//                          static_cast<long>(MILLION) );
// 
//   // Hamming Weight 2, r = 20
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, long > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Weight 2 r=20"} ),
//          HammingWeight2, 30L, 100L * MILLION, 20, 10,
//                          static_cast<long>(MILLION) );
// 
//   // Hamming Correlation, L = 30
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=30"} ),
//          HammingCorr, 1L, 500L * MILLION, 0, 30, 30 );
// 
//   // Hamming Correlation, L = 300
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=300"} ),
//          HammingCorr, 1L, 50L * MILLION, 0, 30, 10 * 30 );
// 
//   // Hamming Correlation, L = 1200
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int> >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Correlation L=1200"} ),
//          HammingCorr, 1L, 10L * MILLION, 0, 30, 40 * 30 );
// 
//   // Hamming independence, L = 30, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=30 r=0"} ),
//          HammingIndep, 1L, 300L * MILLION, 0, 30, 30, 0 );
// 
//   // Hamming independence, L = 30, r = 20
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=30 r=20"} ),
//          HammingIndep, 1L, 100L * MILLION, 20, 10, 30, 0 );
// 
//   // Hamming independence, L = 300, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=300 r=0"} ),
//          HammingIndep, 1L, 30L * MILLION, 0, 30, 10 * 30, 0 );
// 
//   // Hamming independence, L = 300, r = 20
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=300 r=20"} ),
//          HammingIndep, 1L, 10L * MILLION, 20, 10, 10 * 30, 0 );
// 
//   // Hamming independence, L = 1200, r = 0
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=1200 r=0"} ),
//          HammingIndep, 1L, 10L * MILLION, 0, 30, 40 * 30, 0 );
// 
//   // Hamming independence, L = 1200, r = 20
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Hamming Independence L=1200 r=20"} ),
//          HammingIndep, 1L, static_cast<long>(MILLION), 20, 10, 40 * 30, 0 );
// 
//   // String Run, r = 0
//   add< TestU01< sstring_Res3, sstring_CreateRes3, sstring_DeleteRes3,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"String Run NRuns r=0",
//                                          "String Run NBits r=0"} ),
//          StringRun, 1L, 1L * BILLION, 0, 30 );
// 
//   // String Run, r = 20
//   add< TestU01< sstring_Res3, sstring_CreateRes3, sstring_DeleteRes3,
//                 long, long, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"String Run NRuns r=20",
//                                          "String Run NBits r=20"} ),
//          StringRun, 1L, 1L * BILLION, 20, 10 );
// 
//   // Autocorrelation, d = 1, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=1 r=0"} ),
//          AutoCor, 10L, 30L + BILLION, 0, 30, 1 );
// 
//   // Autocorrelation, d = 1, r = 20
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=1 r=20"} ),
//          AutoCor, 5L, 1L + BILLION, 20, 10, 1 );
// 
//   // Autocorrelation, d = 30, r = 0
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=30 r=0"} ),
//          AutoCor, 10L, 31L + BILLION, 0, 30, 30 );
// 
//   // Autocorrelation, d = 10, r = 20
//   add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
//                 long, long, int, int, int > >
//        ( tests, id, gen, rng, StatTest::Names( {"Autocorrelation d=10 r=20"} ),
//          AutoCor, 5L, 11L + BILLION, 20, 10, 10 );
}
