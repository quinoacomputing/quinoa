// *****************************************************************************
/*!
  \file      src/RNGTest/BigCrush.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Class re-creating the TestU01 library's BigCrush battery
  \details   Class re-creating the TestU01 library's BigCrush battery.
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
#include "BigCrush.h"

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::BigCrush;

void
BigCrush::addTests( std::vector< std::function< StatTest() > >& tests,
                    tk::ctr::RNGType rng,
                    CProxy_TestU01Suite& proxy )
// *****************************************************************************
// Add statistical tests to battery
//! \details This function adds, i.e., registers, all statistical tests to the
//!   test stack corresponding to the TestU01 library's BigCrush battery.
//! \param[in] tests Vector of test constructors
//! \param[in] rng RNG ID enum associated with the RNG to be tested
//! \param[in] proxy Charm++ host proxy to which the tests will call back to
//! \author  J. Bakosi
// *****************************************************************************
{
  // Select test stack
  const auto& stack = g_testStack.TestU01;

  // Find RNG
  auto gen = stack.generator( rng );

  static const long THOUSAND = 1000;
  static const long MILLION = THOUSAND * THOUSAND;
  static const long BILLION = THOUSAND * MILLION;

  // Marsaglia Serial Over, r = 0
  stack.add< TestU01< TestU01Props< tag::SerialOver, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Serial Over r=0"},
             1L, static_cast<long>(BILLION), 0, 256L, 3 );

  // Marsaglia Serial Over, r = 22
  stack.add< TestU01< TestU01Props< tag::SerialOver, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Serial Over r=22"},
             1L, static_cast<long>(BILLION), 22, 256L, 3 );

  // Marsaglia Collision Over, t = 2, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=2 r=0"},
             30L, 20L * MILLION, 0, 1024L * 1024L * 2L, 2 );

  // Marsaglia Collision Over, t = 2, r = 9
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=2 r=9"},
             30L, 20L * MILLION, 9, 1024L * 1024L * 2L, 2 );

  // Marsaglia Collision Over, t = 3, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=3 r=0"},
             30L, 20L * MILLION, 0, 1024L * 16L, 3 );

  // Marsaglia Collision Over, t = 3, r = 16
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=3 r=16"},
             30L, 20L * MILLION, 16, 1024L * 16L, 3 );

  // Marsaglia Collision Over, t = 7, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=7 r=0"},
             30L, 20L * MILLION, 0, 64L, 7 );

  // Marsaglia Collision Over, t = 7, r = 24
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=7 r=24"},
             30L, 20L * MILLION, 24, 64L, 7 );

  // Marsaglia Collision Over, t = 14, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=14 r=0"},
             30L, 20L * MILLION, 0, 8L, 14 );

  // Marsaglia Collision Over, t = 14, r = 27
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=14 r=27"},
             30L, 20L * MILLION, 27, 8L, 14 );

  // Marsaglia Collision Over, t = 21, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=21 r=0"},
             30L, 20L * MILLION, 0, 4L, 21 );

  // Marsaglia Collision Over, t = 21, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=21 r=0"},
             30L, 20L * MILLION, 0, 4L, 21 );

  // Marsaglia Collision Over, t = 21, r = 28
  stack.add< TestU01< TestU01Props< tag::CollisionOver, CProxy_TestU01Suite,
                                    smarsa_Res, smarsa_CreateRes,
                                    smarsa_DeleteRes,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Collision Over t=21 r=28"},
             30L, 20L * MILLION, 28, 4L, 21 );

  #ifdef USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  #if LONG_MAX <= 2147483647L
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             250L, 4L * MILLION, 0, 1073741824L, 2, 1 );
  #else // LONG_MAX
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             100L, 10L * MILLION, 0, 2147483648L, 2, 1 );
  #endif // LONG_MAX

  // Marsaglia's Birthday Spacings, t = 3, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=3 r=0"},
             20L, 20L * MILLION, 0, 2097152L, 3, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 14
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=4 r=14"},
             20L, 30L * MILLION, 14, 65536L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=7 r=0"},
             20L, 20L * MILLION, 0, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 7, r = 7
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=7 r=7"},
             20L, 20L * MILLION, 7, 512L, 7, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 14
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=8 r=14"},
             20L, 30L * MILLION, 14, 256L, 8, 1 );

  // Marsaglia's Birthday Spacings, t = 8, r = 22
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=8 r=22"},
             20L, 30L * MILLION, 22, 256L, 8, 1 );

  // Marsaglia's Birthday Spacings, t = 16, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=16 r=0"},
             20L, 30L * MILLION, 0, 16L, 16, 1 );

  // Marsaglia's Birthday Spacings, t = 16, r = 26
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=16 r=26"},
             20L, 30L * MILLION, 26, 16L, 16, 1 );

  #else // USE_LONGLONG

  // Marsaglia's Birthday Spacings, t = 2, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=2 r=0"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 0, 67108864L, 2, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=4 r=0"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 0, 1024L*8L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 4, r = 16
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=4 r=16"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 16, 1024L*8L, 4, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 0
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=0"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 0, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 5
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=5"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 5, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 10
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=10"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 10, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 15
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=15"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 15, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 20
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=20"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 20, 16L, 13, 1 );

  // Marsaglia's Birthday Spacings, t = 13, r = 26
  stack.add< TestU01< TestU01Props< tag::BirthdaySpacings, CProxy_TestU01Suite,
                                    sres_Poisson, sres_CreatePoisson,
                                    sres_DeletePoisson,
                                    long, long, int, long, int, int > > >
           ( proxy, tests, rng, gen, {"Birthday Spacings t=13 r=26"},
             10L*THOUSAND, static_cast<long>(MILLION/10), 26, 16L, 13, 1 );

  #endif // USE_LONGLONG

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
             30L, 6L * MILLION, 0, 3, 0, 30, 1 );

  // Close Pairs, t = 5
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=5",
                                      "Close Pairs mNP t=5",
                                      "Close Pairs mNP1 t=5",
                                      "Close Pairs mNP2 t=5",
                                      "Close Pairs mNJumps t=5",
                                      "Close Pairs mNP2S t=5"},
             20L, 4L * MILLION, 0, 5, 0, 30, 1 );

  // Close Pairs, t = 9
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=9",
                                      "Close Pairs mNP t=9",
                                      "Close Pairs mNP1 t=9",
                                      "Close Pairs mNP2 t=9",
                                      "Close Pairs mNJumps t=9",
                                      "Close Pairs mNP2S t=9"},
             10L, 3L * MILLION, 0, 9, 0, 30, 1 );

  // Close Pairs, t = 16
  stack.add< TestU01< TestU01Props< tag::ClosePairs, CProxy_TestU01Suite,
                                    snpair_Res, snpair_CreateRes,
                                    snpair_DeleteRes,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Close Pairs NP t=16",
                                      "Close Pairs mNP t=16",
                                      "Close Pairs mNP1 t=16",
                                      "Close Pairs mNP2 t=16",
                                      "Close Pairs mNJumps t=16",
                                      "Close Pairs mNP2S t=16"},
             5L, 2L * MILLION, 0, 16, 0, 30, 1 );

  // Knuth's Simple Poker, d = 8, r = 0
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Simplified Poker d=8 r=0"},
               1L, 400L * MILLION, 0, 8, 8 );

  // Knuth's Simple Poker, d = 8, r = 27
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Simplified Poker d=8 r=27"},
               1L, 400L * MILLION, 27, 8, 8 );

  // Knuth's Simple Poker, d = 32, r = 0
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Simplified Poker d=32 r=0"},
               1L, 100L * MILLION, 0, 32, 32 );

  // Knuth's Simple Poker, d = 32, r = 25
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Simplified Poker d=32 r=25"},
               1L, 100L * MILLION, 25, 32, 32 );

  // Knuth's Coupon Collector, d = 8, r = 0
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Coupon Collector d=8 r=0"},
               1L, 200L * MILLION, 0, 8 );

  // Knuth's Coupon Collector, d = 8, r = 10
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Coupon Collector d=8 r=10"},
               1L, 200L * MILLION, 10, 8 );

  // Knuth's Coupon Collector, d = 8, r = 20
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Coupon Collector d=8 r=20"},
               1L, 200L * MILLION, 20, 8 );

  // Knuth's Coupon Collector, d = 8, r = 27
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Coupon Collector d=8 r=27"},
               1L, 200L * MILLION, 27, 8 );

  // Knuth's Gap, r = 0
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
             ( proxy, tests, rng, gen, {"Gap r=0"},
               1L, static_cast<long>(BILLION/2), 0, 0.0, 1.0/16.0 );

  // Knuth's Gap, r = 25
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
             ( proxy, tests, rng, gen, {"Gap r=25"},
               1L, 300L * MILLION, 25, 0.0, 1.0/32.0 );

  // Knuth's Gap, r = 0, n = .5e+9
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
             ( proxy, tests, rng, gen, {"Gap r=0 n=.5B"},
               1L, static_cast<long>(BILLION / 10), 0, 0.0, 1.0/128.0 );

  // Knuth's Gap, r = 20, n = 10e+6
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
             ( proxy, tests, rng, gen, {"Gap r=20 n=10M"},
               1L, 10L * MILLION, 20, 0.0, 1.0/1024.0 );

  // Knuth's Run, r = 0
  stack.add< TestU01< TestU01Props< tag::Run, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Run r=0"},
               5L, static_cast<long>(BILLION), 0, 0 );

  // Knuth's Run, r = 15
  stack.add< TestU01< TestU01Props< tag::Run, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Run r=15"},
               10L, static_cast<long>(BILLION), 15, 1 );

  // Knuth's Permutation, t = 3
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Permutation t=3"},
               1L, static_cast<long>(BILLION), 5, 3 );

  // Knuth's Permutation, t = 5
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Permutation t=5"},
               1L, static_cast<long>(BILLION), 5, 5 );

  // Knuth's Permutation, t = 7
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Permutation t=7"},
               1L, static_cast<long>(BILLION / 2), 5, 7 );

  // Knuth's Permutation, t = 10
  stack.add< TestU01< TestU01Props< tag::Permutation, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Permutation t=10"},
               1L, static_cast<long>(BILLION / 2), 10, 10 );

  // Knuth's Collision with permutations, r = 0
  stack.add< TestU01< TestU01Props< tag::CollisionPermut, CProxy_TestU01Suite,
                                    sknuth_Res2, sknuth_CreateRes2,
                                    sknuth_DeleteRes2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Collision w. Permutations r=0"},
               20L, 20L * MILLION, 0, 14 );

  // Knuth's Collision with permutations, r = 10
  stack.add< TestU01< TestU01Props< tag::CollisionPermut, CProxy_TestU01Suite,
                                    sknuth_Res2, sknuth_CreateRes2,
                                    sknuth_DeleteRes2,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Collision w. Permutations r=10"},
               20L, 20L * MILLION, 10, 14 );

  // Knuth's Maximum-of-t, t = 8
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Maximum-of-t t=8",
                                        "Maximum-of-t Anderson-Darling t=8"},
               40L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 8,
                 gofw_Sum, gofw_AD );

  // Knuth's Maximum-of-t, t = 16
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Maximum-of-t t=16",
                                        "Maximum-of-t Anderson-Darling t=16"},
               30L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 16,
                 gofw_Sum, gofw_AD );

  // Knuth's Maximum-of-t, t = 24
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Maximum-of-t t=24",
                                        "Maximum-of-t Anderson-Darling t=24"},
               20L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 24,
                 gofw_Mean, gofw_Mean );

  // Knuth's Maximum-of-t, t = 32
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Maximum-of-t t=32",
                                        "Maximum-of-t Anderson-Darling t=32"},
               20L, 10L * MILLION, 0, static_cast<int>(MILLION / 10), 32,
                 gofw_Mean, gofw_Mean );

  // Sample Products, t = 8
  stack.add< TestU01< TestU01Props< tag::SampleProd, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Sample Products t=8"},
               4L, 10L * MILLION, 0, 8 );

  // Sample Products, t = 16
  stack.add< TestU01< TestU01Props< tag::SampleProd, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Sample Products t=16"},
               20L, 10L * MILLION, 0, 16 );

  // Sample Products, t = 24
  stack.add< TestU01< TestU01Props< tag::SampleProd, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Sample Products t=24"},
               20L, 10L * MILLION, 0, 24 );

  // Sample Mean, r = 0
  stack.add< TestU01< TestU01Props< tag::SampleMean, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int > > >
             ( proxy, tests, rng, gen, {"Sample Mean r=0"},
               20L * MILLION, 30L, 0 );

  // Sample Mean, r = 10
  stack.add< TestU01< TestU01Props< tag::SampleMean, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int > > >
             ( proxy, tests, rng, gen, {"Sample Mean r=10"},
               20L * MILLION, 30L, 10 );

  // Sample Autorrelation, k = 1
  stack.add< TestU01< TestU01Props< tag::SampleCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Sample Autorrelation k=1"},
               1L, 2L * BILLION, 0, 1 );

  // Sample Autorrelation, k = 2
  stack.add< TestU01< TestU01Props< tag::SampleCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int > > >
             ( proxy, tests, rng, gen, {"Sample Autorrelation k=2"},
               1L, 2L * BILLION, 0, 2 );

  // Maurer's "universal" test, r = 0
  stack.add< TestU01< TestU01Props< tag::AppearanceSpacings, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Appearance Spacings r=0"},
               1L, 10L * MILLION, 400L * MILLION, 0, 30, 15 );

  // Maurer's "universal" test, r = 0
  stack.add< TestU01< TestU01Props< tag::AppearanceSpacings, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Appearance Spacings r=0"},
               1L, 10L * MILLION, static_cast<long>(BILLION), 0,
                             3, 15 );

  // Maurer's "universal" test, r = 27
  stack.add< TestU01< TestU01Props< tag::AppearanceSpacings, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, long, int, int, int > > >
             ( proxy, tests, rng, gen, {"Appearance Spacings r=27"},
               1L, 10L * MILLION, static_cast<long>(BILLION), 27,
                             3, 15 );

  // Weight Distribution, r = 0
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=0"},
               1L, 20L * MILLION, 0, 256L, 0.0, 0.25 );

  // Weight Distribution, r = 20
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=20"},
               1L, 20L * MILLION, 20, 256L, 0.0, 0.25 );

  // Weight Distribution, r = 28
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=28"},
               1L, 20L * MILLION, 28, 256L, 0.0, 0.25 );

  // Weight Distribution, r = 0, beta = 1/16
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=0 beta=1/16"},
               1L, 20L * MILLION, 0, 256L, 0.0, 0.0625 );

  // Weight Distribution, r = 10
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=10"},
               1L, 20L * MILLION, 10, 256L, 0.0, 0.0625 );

  // Weight Distribution, r = 26
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
             ( proxy, tests, rng, gen, {"Weight Distribution r=26"},
               1L, 20L * MILLION, 26, 256L, 0.0, 0.0625 );

  // Sum Collector
  stack.add< TestU01< TestU01Props< tag::SumCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double > > >
             ( proxy, tests, rng, gen, {"Sum Collector"},
               1L, 500L * MILLION, 0, 10.0 );

  // Marsaglia's Matrix Rank, L = 30, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Matrix Rank L=30 r=0"},
               10L, static_cast<long>(MILLION), 0, 5, 30, 30 );

  // Marsaglia's Matrix Rank, L = 30, r = 25
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Matrix Rank L=30 r=25"},
               10L, static_cast<long>(MILLION), 25, 5, 30, 30 );

  // Marsaglia's Matrix Rank, L = 1000, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Matrix Rank L=1000 r=0"},
               1L, 5L * THOUSAND, 0, 4, 1000, 1000 );

  // Marsaglia's Matrix Rank, L = 1000, r = 26
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
             ( proxy, tests, rng, gen, {"Matrix Rank L=1000 r=26"},
               1L, 5L * THOUSAND, 26, 4, 1000, 1000 );

  // Marsaglia's Matrix Rank, L = 5000, r = 15
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank L=5000 r=15"},
             1L, 80L, 15, 15, 5000, 5000 );

  // Marsaglia's Matrix Rank, L = 5000, r = 0
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank L=5000 r=0"},
             1L, 80L, 0, 30, 5000, 5000 );

  // Marsaglia's Modified Savir
  stack.add< TestU01< TestU01Props< tag::Savir2, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, int > > >
           ( proxy, tests, rng, gen, {"Savir2"},
             10L, 10L * MILLION, 10, 1024L * 1024L, 30 );

  // Marsaglia's greatest common divisor
  stack.add< TestU01< TestU01Props< tag::GCD, CProxy_TestU01Suite,
                                    smarsa_Res2, smarsa_CreateRes2,
                                    smarsa_DeleteRes2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"GCD"},
             10L, 50L * MILLION, 0, 30 );

  // Random Walk 1, L = 50, r = 0
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=50 r=0",
                                      "Random Walk 1 Stat M L=50 r=0",
                                      "Random Walk 1 Stat J L=50 r=0",
                                      "Random Walk 1 Stat R L=50 r=0",
                                      "Random Walk 1 Stat C L=50 r=0"},
             1L, 100L * MILLION, 0, 5, 50L, 50L );

  // Random Walk 1, L = 50, r = 25
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=50 r=25",
                                      "Random Walk 1 Stat M L=50 r=25",
                                      "Random Walk 1 Stat J L=50 r=25",
                                      "Random Walk 1 Stat R L=50 r=25",
                                      "Random Walk 1 Stat C L=50 r=25"},
             1L, 100L * MILLION, 25, 5, 50L, 50L );

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
             1L, 10L * MILLION, 0, 10, 1000L, 1000L );

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
             1L, 10L * MILLION, 20, 10, 1000L, 1000L );

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
             1L, static_cast<long>(MILLION), 0, 15, 10000L, 10000L );

  // Random Walk 1, L = 10000, r = 15
  stack.add< TestU01< TestU01Props< tag::RandomWalk1, CProxy_TestU01Suite,
                                    swalk_Res, swalk_CreateRes,
                                    swalk_DeleteRes,
                                    long, long, int, int, long, long > > >
           ( proxy, tests, rng, gen, {"Random Walk 1 Stat H L=10000 r=15",
                                      "Random Walk 1 Stat M L=10000 r=15",
                                      "Random Walk 1 Stat J L=10000 r=15",
                                      "Random Walk 1 Stat R L=10000 r=15",
                                      "Random Walk 1 Stat C L=10000 r=15"},
             1L, static_cast<long>(MILLION), 15, 15, 10000L, 10000L );

  // Linear Complexity, r = 0
  stack.add< TestU01< TestU01Props< tag::LinearComp, CProxy_TestU01Suite,
                                    scomp_Res, scomp_CreateRes,
                                    scomp_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Linear Complexity Jump r=0",
                                      "Linear Complexity Size r=0"},
             1L, 400L * THOUSAND + 20, 0, 1 );

  // Linear Complexity, r = 29
  stack.add< TestU01< TestU01Props< tag::LinearComp, CProxy_TestU01Suite,
                                    scomp_Res, scomp_CreateRes,
                                    scomp_DeleteRes,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Linear Complexity Jump r=29",
                                      "Linear Complexity Size r=29"},
             1L, 400L * THOUSAND + 20, 29, 1 );

  // Lempel-Ziv Compressibility, r = 0
  stack.add< TestU01< TestU01Props< tag::LempelZiv, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Lempel-Ziv Compressibility r=0"},
             10L, 27, 0, 30 );

  // Lempel-Ziv Compressibility, r = 15
  stack.add< TestU01< TestU01Props< tag::LempelZiv, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Lempel-Ziv Compressibility r=15"},
             10L, 27, 15, 15 );

  // Fourier3, r = 0
  stack.add< TestU01< TestU01Props< tag::Fourier3, CProxy_TestU01Suite,
                                    sspectral_Res, sspectral_CreateRes,
                                    sspectral_DeleteRes,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Fourier 3 r=0"},
             100L * THOUSAND, 14, 0, 3 );

  // Fourier3, r = 20
  stack.add< TestU01< TestU01Props< tag::Fourier3, CProxy_TestU01Suite,
                                    sspectral_Res, sspectral_CreateRes,
                                    sspectral_DeleteRes,
                                    long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Fourier 3 r=27"},
             100L * THOUSAND, 14, 27, 3 );

  // Longes Head Run, r = 0
  stack.add< TestU01< TestU01Props< tag::LongestHeadRun, CProxy_TestU01Suite,
                                    sstring_Res2, sstring_CreateRes2,
                                    sstring_DeleteRes2,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Longest Head Run Chi r=0",
                                      "Longest Head Run Disc r=0"},
             1L, 1000L, 0, 3, 20L + 10L * MILLION );

  // Longes Head Run, r = 27
  stack.add< TestU01< TestU01Props< tag::LongestHeadRun, CProxy_TestU01Suite,
                                    sstring_Res2, sstring_CreateRes2,
                                    sstring_DeleteRes2,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Longest Head Run Chi r=27",
                                      "Longest Head Run Disc r=27"},
             1L, 1000L, 27, 3, 20L + 10L * MILLION );

  // Periods In Strings, r = 0
  stack.add< TestU01< TestU01Props< tag::PeriodsInStrings, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Periods In Strings r=0"},
             10L, BILLION / 2, 0, 10 );

  // Periods In Strings, r = 20
  stack.add< TestU01< TestU01Props< tag::PeriodsInStrings, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Periods In Strings r=20"},
             10L, BILLION / 2, 20, 10 );

  // Hamming Weight 2, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingWeight2, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Hamming Weight 2 r=0"},
             10L, static_cast<long>(BILLION), 0, 3,
                         static_cast<long>(MILLION) );

  // Hamming Weight 2, r = 27
  stack.add< TestU01< TestU01Props< tag::HammingWeight2, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, long > > >
           ( proxy, tests, rng, gen, {"Hamming Weight 2 r=27"},
             10L, static_cast<long>(BILLION), 27, 3,
                         static_cast<long>(MILLION) );

  // Hamming Correlation, L = 30
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int> > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=30"},
             1L, static_cast<long>(BILLION), 10, 10, 30 );

  // Hamming Correlation, L = 300
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int> > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=300"},
             1L, 100L * MILLION, 10, 10, 10 * 30 );

  // Hamming Correlation, L = 1200
  stack.add< TestU01< TestU01Props< tag::HammingCorr, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int> > >
           ( proxy, tests, rng, gen, {"Hamming Correlation L=1200"},
             1L, 100L * MILLION, 10, 10, 40 * 30 );

  // Hamming independence, L = 30, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=30 r=0"},
             10L, 30L * MILLION, 0, 3, 30, 0 );

  // Hamming independence, L = 30, r = 27
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=30 r=27"},
             10L, 30L * MILLION, 27, 3, 30, 0 );

  // Hamming independence, L = 300, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=300 r=0"},
             1L, 30L * MILLION, 0, 4, 10 * 30, 0 );

  // Hamming independence, L = 300, r = 26
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=300 r=26"},
             1L, 30L * MILLION, 26, 4, 10 * 30, 0 );

  // Hamming independence, L = 1200, r = 0
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=1200 r=0"},
             1L, 10L * MILLION, 0, 5, 40 * 30, 0 );

  // Hamming independence, L = 1200, r = 25
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence L=1200 r=25"},
             1L, 10L * MILLION, 25, 5, 40 * 30, 0 );

  // String Run, r = 0
  stack.add< TestU01< TestU01Props< tag::StringRun, CProxy_TestU01Suite,
                                    sstring_Res3, sstring_CreateRes3,
                                    sstring_DeleteRes3,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"String Run NRuns r=0",
                                      "String Run NBits r=0"},
             1L, 2L * BILLION, 0, 3 );

  // String Run, r = 27
  stack.add< TestU01< TestU01Props< tag::StringRun, CProxy_TestU01Suite,
                                    sstring_Res3, sstring_CreateRes3,
                                    sstring_DeleteRes3,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"String Run NRuns r=27",
                                      "String Run NBits r=27"},
             1L, 2L * BILLION, 27, 3 );

  // Autocorrelation, d = 1, r = 0
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=1 r=0"},
             10L, 30L + BILLION, 0, 3, 1 );

  // Autocorrelation, d = 3, r = 0
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=3 r=0"},
             10L, 30L + BILLION, 0, 3, 3 );

  // Autocorrelation, d = 1, r =27
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=1 r=27"},
             10L, 30L + BILLION, 27, 3, 1 );

  // Autocorrelation, d = 3, r = 27
  stack.add< TestU01< TestU01Props< tag::AutoCorr, CProxy_TestU01Suite,
                                    sres_Basic, sres_CreateBasic,
                                    sres_DeleteBasic,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Autocorrelation d=3 r=27"},
             10L, 30L + BILLION, 27, 3, 3 );
}
