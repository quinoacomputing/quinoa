//******************************************************************************
/*!
  \file      src/Control/RNGTest/Tags.h
  \author    J. Bakosi
  \date      Fri 04 Jul 2014 06:59:27 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     RNGTest's input deck tags
  \details   RNGTest's input deck tags
*/
//******************************************************************************
#ifndef RNGTestInputDeckTags_h
#define RNGTestInputDeckTags_h

namespace rngtest {
namespace tag {

struct io {};
struct control {};
struct title {};
struct selected {};
struct battery {};
struct generator {};
struct cmd {};
struct param {};

struct BirthdaySpacings {};
struct Collision {};
struct RandomWalk1 {};
struct Gap {};
struct SimplePoker {};
struct CouponCollector {};
struct MaxOft {};
struct WeightDistrib {};
struct MatrixRank {};
struct HammingIndep {};
struct SerialOver {};
struct CollisionOver {};
struct ClosePairs {};
struct ClosePairsBitMatch {};
struct Run {};
struct Permutation {};
struct CollisionPermut {};
struct SampleProd {};
struct SampleMean {};
struct SampleCorr {};
struct AppearanceSpacings {};
struct SumCollector {};
struct Savir2 {};
struct GCD {};
struct LinearComp {};
struct LempelZiv {};
struct Fourier3 {};
struct LongestHeadRun {};
struct PeriodsInStrings {};
struct HammingWeight2 {};
struct HammingCorr {};
struct StringRun {};
struct AutoCorr {};

} // tag::
} // rngtest::

#endif // RNGTestInputDeckTags_h
