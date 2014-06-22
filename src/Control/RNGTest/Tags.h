//******************************************************************************
/*!
  \file      src/Control/RNGTest/Tags.h
  \author    J. Bakosi
  \date      Sat 21 Jun 2014 04:56:14 PM MDT
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

// TestU01 statistical test tags
struct BirthdaySpacings {};
struct Collision {};
struct SerialOver {};
struct RandomWalk1 {};

} // tag::
} // rngtest::

#endif // RNGTestInputDeckTags_h
