// *****************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.C
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
  \brief     Class re-creating the TestU01 library's SmallCrush battery
  \details   Class re-creating the TestU01 library's SmallCrush battery.
*/
// *****************************************************************************

#include "NoWarning/charm.h"
#include "NoWarning/pup_stl.h"
#include "NoWarning/pup.h"

extern "C" {
  #include <gdef.h>
  #include <gofw.h>
  #include <sknuth.h>
  #include <sres.h>
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
#include "SmallCrush.h"

namespace rngtest {

extern TestStack g_testStack;

} // rngtest::

using rngtest::SmallCrush;

void
SmallCrush::addTests( std::vector< std::function< StatTest() > >& tests,
                      tk::ctr::RNGType rng,
                      CProxy_TestU01Suite& proxy )
// *****************************************************************************
// Add statistical tests to battery
//! \details This function adds, i.e., registers, all statistical tests to the
//!   test stack corresponding to the TestU01 library's SmallCrush battery.
//! \param[in] tests Vector of test constructors to add tests to
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

  // Knuth's Gap
  stack.add< TestU01< TestU01Props< tag::Gap, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, double, double > > >
           ( proxy, tests, rng, gen, {"Gap"},
             1L, MILLION / 5, 22, 0.0, 0.00390625 );

  // Knuth's Simple Poker
  stack.add< TestU01< TestU01Props< tag::SimplePoker, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int > > >
           ( proxy, tests, rng, gen, {"Simplified Poker"},
             1L, 2L * MILLION / 5, 24, 64, 64 );

  // Knuth's Coupon Collector
  stack.add< TestU01< TestU01Props< tag::CouponCollector, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int > > >
           ( proxy, tests, rng, gen, {"Coupon Collector"},
             1L, MILLION / 2, 26, 16 );

  // Knuth's Maximum-of-t
  stack.add< TestU01< TestU01Props< tag::MaxOft, CProxy_TestU01Suite,
                                    sknuth_Res1, sknuth_CreateRes1,
                                    sknuth_DeleteRes1,
                                    long, long, int, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Maximum-of-t",
                                      "Maximum-of-t Anderson-Darling"},
             1L, 2L * MILLION, 0, static_cast<int>(MILLION / 10), 6, gofw_Mean,
             gofw_Mean );

  // Weight Distribution
  stack.add< TestU01< TestU01Props< tag::WeightDistrib, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, long, double, double > > >
           ( proxy, tests, rng, gen, {"Weight Distribution"},
             1L, MILLION / 5, 27, 256L, 0.0, 0.125 );

  // Marsaglia's Matrix Rank
  stack.add< TestU01< TestU01Props< tag::MatrixRank, CProxy_TestU01Suite,
                                    sres_Chi2, sres_CreateChi2,
                                    sres_DeleteChi2,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Matrix Rank"},
             1L, 20L * THOUSAND, 20, 10, 60, 60 );

  // Hamming independence
  stack.add< TestU01< TestU01Props< tag::HammingIndep, CProxy_TestU01Suite,
                                    sstring_Res, sstring_CreateRes,
                                    sstring_DeleteRes,
                                    long, long, int, int, int, int > > >
           ( proxy, tests, rng, gen, {"Hamming Independence"},
             1L, MILLION/2, 20, 10, 300, 0 );

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
