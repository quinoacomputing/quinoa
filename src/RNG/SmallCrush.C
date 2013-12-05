//******************************************************************************
/*!
  \file      src/RNG/SmallCrush.C
  \author    J. Bakosi
  \date      Wed 04 Dec 2013 09:24:42 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     SmallCrush battery
  \details   SmallCrush battery
*/
//******************************************************************************

#include <SmallCrush.h>
#include <TestU01.h>

using rngtest::SmallCrush;

SmallCrush::SmallCrush(const Base& base) :
  TestU01Suite( base,
    tk::Option< ctr::Battery >().name( rngtest::ctr::BatteryType::SMALLCRUSH ) )
//******************************************************************************
//  Constructor
//! \author  J. Bakosi
//******************************************************************************
{
  setupRNGs( *this );
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
