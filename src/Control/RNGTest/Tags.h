//******************************************************************************
/*!
  \file      src/Control/RNGTest/Tags.h
  \author    J. Bakosi
  \date      Sat 28 Jun 2014 09:25:29 PM MDT
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

} // tag::
} // rngtest::

#endif // RNGTestInputDeckTags_h
