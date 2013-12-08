//******************************************************************************
/*!
  \file      src/RNG/BigCrush.C
  \author    J. Bakosi
  \date      Sat 07 Dec 2013 08:58:26 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     BigCrush battery
  \details   BigCrush battery
*/
//******************************************************************************

#include <BigCrush.h>
#include <TestU01.h>

using rngtest::BigCrush;

BigCrush::BigCrush(const Base& base) :
  TestU01Suite( base,
    tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::BIGCRUSH ) )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  setupRNGs( *this );
}

void
BigCrush::addTests( const quinoa::ctr::RNGType& rng, const Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
  // Marsaglia Serial Over, r = 0
  add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Serial Over r=0"} ),
       SerialOver, 1L, static_cast<long>(BILLION), 0, 256L, 3 );

  // Marsaglia Serial Over, r = 22
  add< TestU01< sres_Basic, sres_CreateBasic, sres_DeleteBasic,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Serial Over r=22"} ),
       SerialOver, 1L, static_cast<long>(BILLION), 22, 256L, 3 );

  // Marsaglia Collision Over, t = 2, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=2 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 1024L * 1024L * 2L, 2 );

  // Marsaglia Collision Over, t = 2, r = 9
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=2 r=9"} ),
       CollisionOver, 30L, 20L * MILLION, 9, 1024L * 1024L * 2L, 2 );

  // Marsaglia Collision Over, t = 3, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=3 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 1024L * 16L, 3 );

  // Marsaglia Collision Over, t = 3, r = 16
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=3 r=16"} ),
       CollisionOver, 30L, 20L * MILLION, 16, 1024L * 16L, 3 );

  // Marsaglia Collision Over, t = 7, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=7 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 64L, 7 );

  // Marsaglia Collision Over, t = 7, r = 24
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=7 r=24"} ),
       CollisionOver, 30L, 20L * MILLION, 24, 64L, 7 );

  // Marsaglia Collision Over, t = 14, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=14 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 8L, 14 );

  // Marsaglia Collision Over, t = 14, r = 27
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=14 r=27"} ),
       CollisionOver, 30L, 20L * MILLION, 27, 8L, 14 );

  // Marsaglia Collision Over, t = 21, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=21 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 4L, 21 );

  // Marsaglia Collision Over, t = 21, r = 0
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=21 r=0"} ),
       CollisionOver, 30L, 20L * MILLION, 0, 4L, 21 );

  // Marsaglia Collision Over, t = 21, r = 28
  add< TestU01< smarsa_Res, smarsa_CreateRes, smarsa_DeleteRes,
                long, long, int, long, int > >
     ( gen, rng, StatTest::Names( {"Collision Over t=21 r=28"} ),
       CollisionOver, 30L, 20L * MILLION, 28, 4L, 21 );

  #ifdef USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  #if LONG_MAX <= 2147483647L
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 250L, 4L * MILLION, 0, 1073741824L, 2, 1 );
  #else // LONG_MAX
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 100L, 10L * MILLION, 0, 2147483648L, 2, 1 );
  #endif // LONG_MAX

  // Marsaglia's Birthday Spacings, t = 3, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=3 r=0"} ),
       BirthdaySpacings, 20L, 20L * MILLION, 0, 2097152L, 3, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 14
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=14"} ),
       BirthdaySpacings, 20L, 30L * MILLION, 14, 65536L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=0"} ),
       BirthdaySpacings, 20L, 20L * MILLION, 0, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 7
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=7 r=7"} ),
       BirthdaySpacings, 20L, 20L * MILLION, 7, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 14
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=14"} ),
       BirthdaySpacings, 20L, 30L * MILLION, 14, 256L, 8, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 22
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=8 r=22"} ),
       BirthdaySpacings, 20L, 30L * MILLION, 22, 256L, 8, 1 );

  // Marsaglia's Birthday Spacings, t = 16, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=16 r=0"} ),
       BirthdaySpacings, 20L, 30L * MILLION, 0, 16L, 16, 1 );

  // Marsaglia's Birthday Spacings, t = 16, r = 26
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=16 r=26"} ),
       BirthdaySpacings, 20L, 30L * MILLION, 26, 16L, 16, 1 );

  #else // USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=2 r=0"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
                         67108864L, 2, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=0"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
                         1024L * 8L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 16
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=4 r=16"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 16,
                         1024L * 8L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 0
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=0"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 0,
                         16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 5
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=5"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 5,
                         16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 10
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=10"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 10,
                         16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 15
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=15"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 15,
                         16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 20
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=20"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 20,
                         16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 26
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
     ( gen, rng, StatTest::Names( {"Birthday Spacings t=13 r=26"} ),
       BirthdaySpacings, 10L * THOUSAND, static_cast<long>(MILLION / 10), 26,
                         16L, 13, 1 );

  #endif // USE_LONGLONG

  // Close Pairs, t = 3
  add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
                long, long, int, int, int, int, int > >
     ( gen, rng, StatTest::Names( {"Close Pairs NP t=3",
                                   "Close Pairs mNP t=3",
                                   "Close Pairs mNP1 t=3",
                                   "Close Pairs mNP2 t=3",
                                   "Close Pairs mNJumps t=3",
                                   "Close Pairs mNP2S t=3"} ),
       ClosePairs, 30L, 6L * MILLION, 0, 3, 0, 30, 1 );

  // Close Pairs, t = 5
  add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
                long, long, int, int, int, int, int > >
     ( gen, rng, StatTest::Names( {"Close Pairs NP t=5",
                                   "Close Pairs mNP t=5",
                                   "Close Pairs mNP1 t=5",
                                   "Close Pairs mNP2 t=5",
                                   "Close Pairs mNJumps t=5",
                                   "Close Pairs mNP2S t=5"} ),
       ClosePairs, 20L, 4L * MILLION, 0, 5, 0, 30, 1 );

  // Close Pairs, t = 9
  add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
                long, long, int, int, int, int, int > >
     ( gen, rng, StatTest::Names( {"Close Pairs NP t=9",
                                   "Close Pairs mNP t=9",
                                   "Close Pairs mNP1 t=9",
                                   "Close Pairs mNP2 t=9",
                                   "Close Pairs mNJumps t=9",
                                   "Close Pairs mNP2S t=9"} ),
       ClosePairs, 10L, 3L * MILLION, 0, 9, 0, 30, 1 );

  // Close Pairs, t = 16
  add< TestU01< snpair_Res, snpair_CreateRes, snpair_DeleteRes,
                long, long, int, int, int, int, int > >
     ( gen, rng, StatTest::Names( {"Close Pairs NP t=16",
                                   "Close Pairs mNP t=16",
                                   "Close Pairs mNP1 t=16",
                                   "Close Pairs mNP2 t=16",
                                   "Close Pairs mNJumps t=16",
                                   "Close Pairs mNP2S t=16"} ),
       ClosePairs, 5L, 2L * MILLION, 0, 16, 0, 30, 1 );

  // Knuth's Simple Poker, d = 8, r = 0
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int, int > >
       ( gen, rng, StatTest::Names( {"Simplified Poker d=8 r=0"} ),
         SimpPoker, 1L, 400L * MILLION, 0, 8, 8 );

  // Knuth's Simple Poker, d = 8, r = 27
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int, int > >
       ( gen, rng, StatTest::Names( {"Simplified Poker d=8 r=27"} ),
         SimpPoker, 1L, 400L * MILLION, 27, 8, 8 );

  // Knuth's Simple Poker, d = 32, r = 0
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int, int > >
       ( gen, rng, StatTest::Names( {"Simplified Poker d=32 r=0"} ),
         SimpPoker, 1L, 100L * MILLION, 0, 32, 32 );

  // Knuth's Simple Poker, d = 32, r = 25
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int, int > >
       ( gen, rng, StatTest::Names( {"Simplified Poker d=32 r=25"} ),
         SimpPoker, 1L, 100L * MILLION, 25, 32, 32 );

  // Knuth's Coupon Collector, d = 8, r = 0
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Coupon Collector d=8 r=0"} ),
         CouponCollector, 1L, 200L * MILLION, 0, 8 );

  // Knuth's Coupon Collector, d = 8, r = 10
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Coupon Collector d=8 r=10"} ),
         CouponCollector, 1L, 200L * MILLION, 10, 8 );

  // Knuth's Coupon Collector, d = 8, r = 20
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Coupon Collector d=8 r=20"} ),
         CouponCollector, 1L, 200L * MILLION, 20, 8 );

  // Knuth's Coupon Collector, d = 8, r = 27
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Coupon Collector d=8 r=27"} ),
         CouponCollector, 1L, 200L * MILLION, 27, 8 );

  // Knuth's Gap, r = 0
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, double, double > >
       ( gen, rng, StatTest::Names( {"Gap r=0"} ),
         Gap, 1L, static_cast<long>(BILLION/2), 0, 0.0, 1.0/16.0 );

  // Knuth's Gap, r = 25
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, double, double > >
       ( gen, rng, StatTest::Names( {"Gap r=25"} ),
         Gap, 1L, 300L * MILLION, 25, 0.0, 1.0/32.0 );

  // Knuth's Gap, r = 0, n = .5e+9
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, double, double > >
       ( gen, rng, StatTest::Names( {"Gap r=0 n=.5B"} ),
         Gap, 1L, static_cast<long>(BILLION / 10), 0, 0.0, 1.0/128.0 );

  // Knuth's Gap, r = 20, n = 10e+6
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, double, double > >
       ( gen, rng, StatTest::Names( {"Gap r=20 n=10M"} ),
         Gap, 1L, 10L * MILLION, 20, 0.0, 1.0/1024.0 );

  // Knuth's Run, r = 0
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Run r=0"} ),
         Run, 5L, static_cast<long>(BILLION), 0, 0 );

  // Knuth's Run, r = 15
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Run r=15"} ),
         Run, 10L, static_cast<long>(BILLION), 15, 1 );

  // Knuth's Permutation, t = 3
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Permutation t=3"} ),
         Permutation, 1L, static_cast<long>(BILLION), 5, 3 );

  // Knuth's Permutation, t = 5
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Permutation t=5"} ),
         Permutation, 1L, static_cast<long>(BILLION), 5, 5 );

  // Knuth's Permutation, t = 7
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Permutation t=7"} ),
         Permutation, 1L, static_cast<long>(BILLION / 2), 5, 7 );

  // Knuth's Permutation, t = 10
  add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                long, long, int, int > >
       ( gen, rng, StatTest::Names( {"Permutation t=10"} ),
         Permutation, 1L, static_cast<long>(BILLION / 2), 10, 10 );


}
