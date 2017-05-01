// *****************************************************************************
/*!
  \file      src/RNGTest/Crush.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Class re-creating the TestU01 library's Crush battery
  \details   Class re-creating the TestU01 library's Crush battery.
*/
// *****************************************************************************

#include "NoWarning/charm.h"
#include "NoWarning/pup_stl.h"
#include "NoWarning/pup.h"

extern "C" {
  #include <gdef.h>
  #include <gofw.h>
  #include <scomp.h>
  #include <sknuth.h>
  #include <smarsa.h>
  #include <snpair.h>
  #include <sres.h>
  #include <sspectral.h>
  #include <sstring.h>
  #include <swalk.h>
}

#include "Tags.h"
#include "PUPUtil.h"
#include "TestU01.h"
#include "StatTest.h"
#include "TestStack.h"
#include "TestU01Props.h"
#include "TestU01Stack.h"
#include "Crush.h"

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::Crush;

void
Crush::addTests( std::vector< std::function< StatTest() > >& tests,
                 tk::ctr::RNGType rng,
                 CProxy_TestU01Suite& proxy )
// *****************************************************************************
// Add statistical tests to battery
//! \details This function adds, i.e., registers, all statistical tests to the
//!   test stack corresponding to the TestU01 library's Crush battery.
//! \param[in] tests Vector of test constructors
//! \param[in] rng RNG ID enum associated with the RNG to be tested
//! \param[in] proxy Charm++ host proxy to which the tests will call back to
//! \author  J. Bakosi
// *****************************************************************************
{
  // Select test stack
  const auto& stack = g_testStack.TestU01;

  // Find RNG
  auto gen = stack.generator(rng);

  static const long THOUSAND = 1000;
  static const long MILLION = THOUSAND * THOUSAND;
  static const long BILLION = THOUSAND * MILLION;

  // Marsaglia's Serial Over, t = 2
  stack.add< TestU01< TestU01Props< tag::SerialOver, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Serial Over t=2"},
             1L, 500L * MILLION, 0, 4096L, 2 );

  // Marsaglia's Serial Over, t = 4
  stack.add< TestU01< TestU01Props< tag::SerialOver, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Serial Over t=4"},
             1L, 300L * MILLION, 0, 64L, 4 );

  // Marsaglia's Collision Over, t = 2, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=2 r=0"},
             10L, 10L * MILLION, 0, 1024L * 1024, 2 );

  // Marsaglia's Collision Over, t = 2, r = 10
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=2 r=10"},
             10L, 10L * MILLION, 10, 1024L * 1024, 2 );

  // Marsaglia's Collision Over, t = 4, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=4 r=0"},
             10L, 10L * MILLION, 0, 1024L, 4 );

  // Marsaglia's Collision Over, t = 4, r = 20
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=4 r=20"},
             10L, 10L * MILLION, 20, 1024L, 4 );

  // Marsaglia's Collision Over, t = 8, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=8 r=0"},
             10L, 10L * MILLION, 0, 32L, 8 );

  // Marsaglia's Collision Over, t = 8, r = 25
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=8 r=25"},
             10L, 10L * MILLION, 25, 32L, 8 );

  // Marsaglia's Collision Over, t = 20, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=20 r=0"},
              10L, 10L * MILLION, 0, 4L, 20 );

  // Marsaglia's Collision Over, t = 20, r = 28
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=20 r=28"},
             10L, 10L * MILLION, 28, 4L, 20 );

  #ifdef USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  #if LONG_MAX <= 2147483647L
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             10L, 10L * MILLION, 0, 1073741824L, 2, 1 );
  #else // LONG_MAX
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             5L, 20L * MILLION, 0, 2L*1073741824L, 2, 1 );
  #endif // LONG_MAX

  // Marsaglia's Birthday Spacings, t = 3, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=3 r=0"},
             5L, 20L * MILLION, 0, 2097152L, 3, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=4 r=0"},
             5L, 20L * MILLION, 0, 65536L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=7 r=0"},
             3L, 20L * MILLION, 0, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 7
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=7 r=7"},
             3L, 20L * MILLION, 7, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 14
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=8 r=14"},
             3L, 20L * MILLION, 14, 256L, 8, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 22
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=8 r=22"},
             3L, 20L * MILLION, 22, 256L, 8, 1 );

  #else // USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             200L, 4L * MILLION / 10, 0, 67108864L, 2, 1 );

  // Marsaglia's Birthday Spacings, t = 3, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=3 r=0"},
             100L, 4L * MILLION / 10, 0, 131072L, 3, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=4 r=0"},
             200L, 4L * MILLION / 10, 0, 1024L * 8, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=0"},
             100L, 4L * MILLION / 10, 0, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 10
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=10"},
             100L, 4L * MILLION / 10, 10, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 20
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=20"},
             100L, 4L * MILLION / 10, 20, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 26
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, rng, {"Birthday Spacings t=13 r=26"},
             100L, 4L * MILLION / 10, 26, 16L, 13, 1 );

  #endif // USE_LONGLONG

  // Close Pairs, t = 2
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=2",
                                      "Close Pairs mNP t=2",
                                      "Close Pairs mNP1 t=2",
                                      "Close Pairs mNP2 t=2",
                                      "Close Pairs mNJumps t=2"},
             10L, 2L * MILLION, 0, 2, 0, 30, 0 );

  // Close Pairs, t = 3
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=3",
                                "Close Pairs mNP t=3",
                                "Close Pairs mNP1 t=3",
                                "Close Pairs mNP2 t=3",
                                "Close Pairs mNJumps t=3",
                                "Close Pairs mNP2S t=3"},
             10L, 2L * MILLION, 0, 3, 0, 30, 1 );

  // Close Pairs, t = 7
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=3",
                                      "Close Pairs mNP t=3",
                                      "Close Pairs mNP1 t=3",
                                      "Close Pairs mNP2 t=3",
                                      "Close Pairs mNJumps t=3",
                                      "Close Pairs mNP2S t=3"},
             5L, 2L * MILLION, 0, 7, 0, 30, 1 );

  // Close Pairs Bit Match, t = 2
  stack.add< TestU01< TestU01Props< tag::ClosePairsBitMatch, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs Bit Match t=2"},
             4L, 4L * MILLION, 0, 2 );

  // Close Pairs Bit Match, t = 4
  stack.add< TestU01< TestU01Props< tag::ClosePairsBitMatch, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs Bit Match t=4"},
             2L, 4L * MILLION, 0, 4 );

  // Knuth's Simple Poker, d = 16, r = 0
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Simplified Poker d=16 r=0"},
             1L, 40L * MILLION, 0, 16, 16 );

  // Knuth's Simple Poker, d = 16, r = 26
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Simplified Poker d=16 r=26"},
             1L, 40L * MILLION, 26, 16, 16 );

  // Knuth's Simple Poker, d = 64, r = 0
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Simplified Poker d=64 r=0"},
             1L, 10L * MILLION, 0, 64, 64 );

  // Knuth's Simple Poker, d = 64, r = 24
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Simplified Poker d=64 r=24"},
             1L, 10L * MILLION, 24, 64, 64 );

  // Knuth's Coupon Collector, d = 4, r = 0
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
       ( proxy, tests, rng, gen, {"Coupon Collector d=4 r=0"},
         1L, 40L * MILLION, 0, 4 );

  // Knuth's Coupon Collector, d = 4, r = 28
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Coupon Collector d=4 r=28"},
             1L, 40L * MILLION, 28, 4 );

  // Knuth's Coupon Collector, d = 16, r = 0
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Coupon Collector d=16 r=0"},
             1L, 10L * MILLION, 0, 16 );

  // Knuth's Coupon Collector, d = 16, r = 26
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Coupon Collector d=16 r=26"},
             1L, 10L * MILLION, 26, 16 );

  // Knuth's Gap, r = 0
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
           ( proxy, tests, rng, gen, {"Gap r=0"},
             1L, 100L * MILLION, 0, 0.0, 0.125 );

  // Knuth's Gap, r = 27
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
           ( proxy, tests, rng, gen, {"Gap r=27"},
             1L, 100L * MILLION, 27, 0.0, 0.125 );

  // Knuth's Gap, r = 0, n = 5e+6
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
           ( proxy, tests, rng, gen, {"Gap r=0 n=5M"},
             1L, 5L * MILLION, 0, 0.0, 1.0/256.0 );

  // Knuth's Gap, r = 22, n = 5e+6
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
           ( proxy, tests, rng, gen, {"Gap r=22 n=5M"},
             1L, 5L * MILLION, 22, 0.0, 1.0/256.0 );

  // Knuth's Run, r = 0
  stack.add< TestU01< TestU01Props< tag::Run, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Run r=0"},
             1L, 500L * MILLION, 0, 1 );

  // Knuth's Run, r = 15
  stack.add< TestU01< TestU01Props< tag::Run, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Run r=15"},
             1L, 500L * MILLION, 15, 0 );

  // Knuth's Permutation, r = 0
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Permutation r=0"},
             1L, 50L * MILLION, 0, 10 );

  // Knuth's Permutation, r = 15
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Permutation r=15"},
             1L, 50L * MILLION, 15, 10 );

  // Knuth's Collision with permutations, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionPermut, CProxy_TestU01Suite,
                                    sknuth_Res2, sknuth_CreateRes2,
                                    sknuth_DeleteRes2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Collision w. Permutations r=0"},
             5L, 10L * MILLION, 0, 13 );

  // Knuth's Collision with permutations, r = 15
  stack.add< TestU01< TestU01Props< tag::CollisionPermut, CProxy_TestU01Suite,
                                    sknuth_Res2, sknuth_CreateRes2,
                                    sknuth_DeleteRes2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Collision w. Permutations r=15"},
             5L, 10L * MILLION, 15, 13 );

  // Knuth's Maximum-of-t, t = 5
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Maximum-of-t t=5",
                                      "Maximum-of-t Anderson-Darling t=5"},
             10L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 5, gofw_Sum,
             gofw_AD );

  // Knuth's Maximum-of-t, t = 10
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Maximum-of-t t=10",
                                      "Maximum-of-t Anderson-Darling t=10"},
             5L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 10, gofw_Sum,
             gofw_AD );

  // Knuth's Maximum-of-t, t = 20
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Maximum-of-t t=20",
                                      "Maximum-of-t Anderson-Darling t=20"},
             1L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 20,
             gofw_Mean, gofw_Mean );

  // Knuth's Maximum-of-t, t = 30
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Maximum-of-t t=30",
                                      "Maximum-of-t Anderson-Darling t=30"},
             1L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 30,
             gofw_Mean, gofw_Mean );

  // Sample Products, t = 10
  stack.add< TestU01< TestU01Props< tag::SampleProd, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Sample Products t=10"},
             1L, 10L * MILLION, 0, 10 );

  // Sample Products, t = 30
  stack. add< TestU01< TestU01Props< tag::SampleProd, CProxy_TestU01Suite,
                                     sres_Basic, sres_CreateBasic,
                                     sres_DeleteBasic,
                                     long, long, int, int > > >
            ( proxy, tests, rng, gen, {"Sample Products t=30"},
              1L, 10L * MILLION, 0, 30 );

  // Sample Mean
  stack.add< TestU01< TestU01Props< tag::SampleMean, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int > > >
           ( proxy, tests, rng, gen, {"Sample Mean"},
             10L * MILLION, 20L, 0 );

  // Sample Autorrelation
  stack.add< TestU01< TestU01Props< tag::SampleCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Sample Autorrelation"},
             1L, 500L * MILLION, 0, 1 );

  // Maurer's "universal" test, r = 0
  stack.add< TestU01< TestU01Props< tag::AppearanceSpacings, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Appearance Spacings r=0"},
             1L, 10L * MILLION, 400L * MILLION, 0, 30, 15 );

  // Maurer's "universal" test, r = 20
  stack.add< TestU01< TestU01Props< tag::AppearanceSpacings, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Appearance Spacings r=20"},
             1L, 10L * MILLION, 100L * MILLION, 20, 10, 15 );

  // Weight Distribution, r = 0
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
           ( proxy, tests, rng, gen, {"Weight Distribution r=0"},
             1L, 2L * MILLION, 0, 256L, 0.0, 0.125 );

  // Weight Distribution, r = 8
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
           ( proxy, tests, rng, gen, {"Weight Distribution r=8"},
             1L, 2L * MILLION, 8, 256L, 0.0, 0.125 );

  // Weight Distribution, r = 16
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
           ( proxy, tests, rng, gen, {"Weight Distribution r=16"},
             1L, 2L * MILLION, 16, 256L, 0.0, 0.125 );

  // Weight Distribution, r = 24
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
           ( proxy, tests, rng, gen, {"Weight Distribution r=24"},
             1L, 2L * MILLION, 24, 256L, 0.0, 0.125 );

  // Sum Collector
  stack.add< TestU01< TestU01Props< tag::SumCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double > > >
           ( proxy, tests, rng, gen, {"Sum Collector"},
             1L, 20L * MILLION, 0, 10.0 );

  // Marsaglia's Matrix Rank, 60 x 60, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank 60x60 r=0"},
             1L, static_cast<long>(MILLION), 0, 30, 2 * 30, 2 * 30 );

  // Marsaglia's Matrix Rank, 60 x 60, r = 20
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank 60x60 r=20"},
             1L, static_cast<long>(MILLION), 20, 10, 2 * 30, 2 * 30 );

  // Marsaglia's Matrix Rank, 300 x 300, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank 300x300 r=0"},
             1L, 50L * THOUSAND, 0, 30, 10 * 30, 10 * 30 );

  // Marsaglia's Matrix Rank, 300 x 300, r = 20
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank 300x300 r=20"},
             1L, 50L * THOUSAND, 20, 10, 10 * 30, 10 * 30 );

  // Marsaglia's Matrix Rank, 1200 x 1200, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
         ( proxy, tests, rng, gen, {"Matrix Rank 1200x1200 r=0"},
           1L, 2L * THOUSAND, 0, 30, 40 * 30, 40 * 30 );

  // Marsaglia's Matrix Rank, 1200 x 1200, r = 20
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank 1200x1200 r=20"},
             1L, 2L * THOUSAND, 20, 10, 40 * 30, 40 * 30 );

  // Marsaglia's Modified Savir
  stack.add< TestU01< TestU01Props< tag::Savir2, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Savir2"},
             1L, 20L * MILLION, 0, 1024L * 1024L, 30 );

  // Marsaglia's greatest common divisor, r = 0
  stack.add< TestU01< TestU01Props< tag::GCD, CProxy_TestU01Suite,
                                    smarsa_Res2, smarsa_CreateRes2,
                                    smarsa_DeleteRes2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"GCD r=0"},
             1L, 100L * MILLION, 0, 30 );

  // Marsaglia's greatest common divisor, r = 10
  stack.add< TestU01< TestU01Props< tag::GCD, CProxy_TestU01Suite,
                                    smarsa_Res2, smarsa_CreateRes2,
                                    smarsa_DeleteRes2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"GCD r=10"},
             1L, 40L * MILLION, 10, 20 );

  // Random Walk 1, L = 90, r = 0
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=90 r=0",
                                      "Random Walk 1 Stat M L=90 r=0",
                                      "Random Walk 1 Stat J L=90 r=0",
                                      "Random Walk 1 Stat R L=90 r=0",
                                      "Random Walk 1 Stat C L=90 r=0"},
             1L, 50L * MILLION, 0, 30, 90L, 90L );

  // Random Walk 1, L = 90, r = 0
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=90 r=20",
                                      "Random Walk 1 Stat M L=90 r=20",
                                      "Random Walk 1 Stat J L=90 r=20",
                                      "Random Walk 1 Stat R L=90 r=20",
                                      "Random Walk 1 Stat C L=90 r=20"},
             1L, 10L * MILLION, 20, 10, 90L, 90L );

  // Random Walk 1, L = 1000, r = 0
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=1000 r=0",
                                      "Random Walk 1 Stat M L=1000 r=0",
                                      "Random Walk 1 Stat J L=1000 r=0",
                                      "Random Walk 1 Stat R L=1000 r=0",
                                      "Random Walk 1 Stat C L=1000 r=0"},
             1L, 5L * MILLION, 0, 30, 1000L, 1000L );

  // Random Walk 1, L = 1000, r = 20
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=1000 r=20",
                                      "Random Walk 1 Stat M L=1000 r=20",
                                      "Random Walk 1 Stat J L=1000 r=20",
                                      "Random Walk 1 Stat R L=1000 r=20",
                                      "Random Walk 1 Stat C L=1000 r=20"},
             1L, static_cast<long>(MILLION), 20, 10, 1000L, 1000L );

  // Random Walk 1, L = 10000, r = 0
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=10000 r=0",
                                      "Random Walk 1 Stat M L=10000 r=0",
                                      "Random Walk 1 Stat J L=10000 r=0",
                                      "Random Walk 1 Stat R L=10000 r=0",
                                      "Random Walk 1 Stat C L=10000 r=0"},
             1L, MILLION / 2, 0, 30, 10000L, 10000L );

  // Random Walk 1, L = 10000, r = 20
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=10000 r=20",
                                      "Random Walk 1 Stat M L=10000 r=20",
                                      "Random Walk 1 Stat J L=10000 r=20",
                                      "Random Walk 1 Stat R L=10000 r=20",
                                      "Random Walk 1 Stat C L=10000 r=20"},
             1L, MILLION / 10, 20, 10, 10000L, 10000L );

  // Linear Complexity, r = 0
  stack.add< TestU01< TestU01Props< tag::LinearComp, CProxy_TestU01Suite,
                                    scomp_Res, scomp_CreateRes,
                                    scomp_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Linear Complexity Jump r=0",
                                      "Linear Complexity Size r=0"},
             1L, 120L * THOUSAND, 0, 1 );

  // Linear Complexity, r = 29
  stack.add< TestU01< TestU01Props< tag::LinearComp, CProxy_TestU01Suite,
                                    scomp_Res, scomp_CreateRes,
                                    scomp_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Linear Complexity Jump r=29",
                                      "Linear Complexity Size r=29"},
             1L, 120L * THOUSAND, 29, 1 );

  // Lempel-Ziv Compressibility
  stack.add< TestU01< TestU01Props< tag::LempelZiv, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Lempel-Ziv Compressibility"},
             10L, 25, 0, 30 );

  // Fourier3, r = 0
  stack.add< TestU01< TestU01Props< tag::Fourier3, CProxy_TestU01Suite,
                                    sspectral_Res, sspectral_CreateRes,
                                    sspectral_DeleteRes,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Fourier 3 r=0"},
             50L * THOUSAND, 14, 0, 30 );

  // Fourier3, r = 20
  stack.add< TestU01< TestU01Props< tag::Fourier3, CProxy_TestU01Suite,
                              sspectral_Res, sspectral_CreateRes,
                              sspectral_DeleteRes,
                              long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Fourier 3 r=20"},
             50L * THOUSAND, 14, 20, 10 );

  // Longest Head Run, r = 0
  stack.add< TestU01< TestU01Props< tag::LongestHeadRun, CProxy_TestU01Suite,
                                    sstring_Res2, sstring_CreateRes2,
                                    sstring_DeleteRes2,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Longest Head Run Chi r=0",
                                      "Longest Head Run Disc r=0"},
             1L, 1000L, 0, 30, 20L + 10L * MILLION );

  // Longest Head Run, r = 20
  stack.add< TestU01< TestU01Props< tag::LongestHeadRun, CProxy_TestU01Suite,
                                    sstring_Res2, sstring_CreateRes2,
                                    sstring_DeleteRes2,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Longest Head Run Chi r=20",
                                      "Longest Head Run Disc r=20"},
             1L, 300L, 20, 10, 20L + 10L * MILLION );

  // Periods In Strings, r = 0
  stack.add< TestU01< TestU01Props< tag::PeriodsInStrings, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Periods In Strings r=0"},
             1L, 300L * MILLION, 0, 30 );

  // Periods In Strings, r = 15
  stack.add< TestU01< TestU01Props< tag::PeriodsInStrings, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Periods In Strings r=15"},
             1L, 300L * MILLION, 15, 15 );

  // Hamming Weight 2, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingWeight2, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Hamming Weight 2 r=0"},
             100L, 100L * MILLION, 0, 30, static_cast<long>(MILLION) );

  // Hamming Weight 2, r = 20
  stack.add< TestU01< TestU01Props< tag::HammingWeight2, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Hamming Weight 2 r=20"},
             30L, 100L * MILLION, 20, 10, static_cast<long>(MILLION) );

  // Hamming Correlation, L = 30
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=30"},
             1L, 500L * MILLION, 0, 30, 30 );

  // Hamming Correlation, L = 300
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=300"},
             1L, 50L * MILLION, 0, 30, 10 * 30 );

  // Hamming Correlation, L = 1200
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=1200"},
             1L, 10L * MILLION, 0, 30, 40 * 30 );

  // Hamming independence, L = 30, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=30 r=0"},
             1L, 300L * MILLION, 0, 30, 30, 0 );

  // Hamming independence, L = 30, r = 20
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=30 r=20"},
             1L, 100L * MILLION, 20, 10, 30, 0 );

  // Hamming independence, L = 300, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=300 r=0"},
             1L, 30L * MILLION, 0, 30, 10 * 30, 0 );

  // Hamming independence, L = 300, r = 20
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=300 r=20"},
             1L, 10L * MILLION, 20, 10, 10 * 30, 0 );

  // Hamming independence, L = 1200, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=1200 r=0"},
             1L, 10L * MILLION, 0, 30, 40 * 30, 0 );

  // Hamming independence, L = 1200, r = 20
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=1200 r=20"},
             1L, static_cast<long>(MILLION), 20, 10, 40 * 30, 0 );

  // String Run, r = 0
  stack.add< TestU01< TestU01Props< tag::StringRun, CProxy_TestU01Suite,
                                    sstring_Res3, sstring_CreateRes3,
                                    sstring_DeleteRes3,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"String Run NRuns r=0",
                                      "String Run NBits r=0"},
             1L, 1L * BILLION, 0, 30 );

  // String Run, r = 20
  stack.add< TestU01< TestU01Props< tag::StringRun, CProxy_TestU01Suite,
                                    sstring_Res3, sstring_CreateRes3,
                                    sstring_DeleteRes3,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"String Run NRuns r=20",
                                      "String Run NBits r=20"},
             1L, 1L * BILLION, 20, 10 );

  // Autocorrelation, d = 1, r = 0
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=1 r=0"},
             10L, 30L + BILLION, 0, 30, 1 );

  // Autocorrelation, d = 1, r = 20
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=1 r=20"},
             5L, 1L + BILLION, 20, 10, 1 );

  // Autocorrelation, d = 30, r = 0
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=30 r=0"},
             10L, 31L + BILLION, 0, 30, 30 );

  // Autocorrelation, d = 10, r = 20
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=10 r=20"},
             5L, 11L + BILLION, 20, 10, 10 );
}
