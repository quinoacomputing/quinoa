//******************************************************************************
/*!
  \file      src/RNGTest/SmallCrush.C
  \author    J. Bakosi
  \date      Sat 10 May 2014 10:19:18 AM MDT
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
  assignTests( *this );
}

void
SmallCrush::addTests( std::size_t id,
                      tk::ctr::RNGType& rng,
                      Gen01Ptr& gen )
//******************************************************************************
// Add statistical tests to battery
//! \author  J. Bakosi
//******************************************************************************
{
  // Marsaglia's BirthdaySpacings
  #ifdef USE_LONGLONG
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
       ( id, gen, rng, {"Birthday Spacings"},
         BirthdaySpacings, 1L, 5L * MILLION, 0, 1073741824L, 2, 1 );
  #else
  add< TestU01< sres_Poisson, sres_CreatePoisson, sres_DeletePoisson,
                long, long, int, long, int, int > >
       ( id, gen, rng, {"Birthday Spacings"},
         BirthdaySpacings, 10L, MILLION / 2, 0, 67108864L, 2, 1 );
  #endif
// 
//   // Knuth's Collision
//   add< TestU01< sknuth_Res2, sknuth_CreateRes2, sknuth_DeleteRes2,
//                 long, long, int, long, int > >
//        ( id, gen, rng, StatTest::Names( {"Collision"} ),
//          Collision, 1L, 5L * MILLION, 0, 65536L, 2 );
// 
//   // Knuth's Gap
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, double, double > >
//        ( id, gen, rng, StatTest::Names( {"Gap"} ),
//          Gap, 1L, MILLION / 5, 22, 0.0, 0.00390625 );
// 
//   // Knuth's Simple Poker
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int > >
//        ( id, gen, rng, StatTest::Names( {"Simplified Poker"} ),
//          SimpPoker, 1L, 2L * MILLION / 5, 24, 64, 64 );
// 
//   // Knuth's Coupon Collector
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int > >
//        ( id, gen, rng, StatTest::Names( {"Coupon Collector"} ),
//          CouponCollector, 1L, MILLION / 2, 26, 16 );
// 
//   // Knuth's Maximum-of-t
//   add< TestU01< sknuth_Res1, sknuth_CreateRes1, sknuth_DeleteRes1,
//                 long, long, int, int, int, gofw_TestType, gofw_TestType > >
//        ( id, gen, rng, StatTest::Names( {"Maximum-of-t",
//                                          "Maximum-of-t Anderson-Darling"} ),
//          MaxOft, 1L, 2L * MILLION, 0, static_cast<int>(MILLION / 10), 6,
//                  gofw_Mean, gofw_Mean );
// 
//   // Weight Distribution
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, long, double, double > >
//        ( id, gen, rng, StatTest::Names( {"Weight Distribution"} ),
//          WeightDistrib, 1L, MILLION / 5, 27, 256L, 0.0, 0.125 );
// 
//   // Marsaglia's Matrix Rank
//   add< TestU01< sres_Chi2, sres_CreateChi2, sres_DeleteChi2,
//                 long, long, int, int, int, int > >
//        ( id, gen, rng, StatTest::Names( {"Matrix Rank"} ),
//          MatrixRank, 1L, 20L * THOUSAND, 20, 10, 60, 60 );
// 
//   // Hamming independence
//   add< TestU01< sstring_Res, sstring_CreateRes, sstring_DeleteRes,
//                 long, long, int, int, int, int > >
//        ( id, gen, rng, StatTest::Names( {"Hamming Independence"} ),
//          HammingIndep, 1L, MILLION/2, 20, 10, 300, 0 );
// 
//   // Random Walk 1
//   add< TestU01< swalk_Res, swalk_CreateRes, swalk_DeleteRes,
//                 long, long, int, int, long, long > >
//        ( id, gen, rng, StatTest::Names( {"Random Walk 1 Stat H",
//                                          "Random Walk 1 Stat M",
//                                          "Random Walk 1 Stat J",
//                                          "Random Walk 1 Stat R",
//                                          "Random Walk 1 Stat C"} ),
//          RandomWalk1, 1L, static_cast<long>(MILLION), 0, 30, 150L, 150L );
}
