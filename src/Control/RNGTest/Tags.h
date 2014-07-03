//******************************************************************************
/*!
  \file      src/Control/RNGTest/Tags.h
  \author    J. Bakosi
  \date      Thu 03 Jul 2014 07:59:06 AM MDT
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
struct SerialOver {};
struct RandomWalk1 {};
struct Gap {};
struct SimplePoker {};
struct CouponCollector {};
struct MaxOft {};
struct WeightDistrib {};
struct MatrixRank {};
struct HammingIndep {};

} // tag::
} // rngtest::

#endif // RNGTestInputDeckTags_h
