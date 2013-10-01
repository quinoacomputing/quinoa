//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Tags.h
  \author    J. Bakosi
  \date      Mon 30 Sep 2013 09:39:00 PM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck tags
  \details   All control tags used to build a nested tagged tagged tuple
*/
//******************************************************************************
#ifndef QuinoaInputDeckTags_h
#define QuinoaInputDeckTags_h

namespace quinoa {
namespace ctr {

struct geometry {};
struct physics {};
struct position {};
struct mass {};
struct hydro {};
struct energy {};
struct mix {};
struct frequency {};
struct mixrate {};

struct nstep {};
struct term {};
struct dt {};

struct nposition {};
struct ndensity {};
struct nvelocity {};
struct nscalar {};
struct nfrequency {};
struct npar {};

struct tty {};
struct dump {};
struct plot {};
struct pdf {};
struct glob {};

struct control {};
struct input {};
struct geomoutput {};
struct physoutput {};
struct stats {};

struct atwood {};
struct b {};
struct S {};
struct kappa {};
struct c {};
struct c0 {};
struct c1 {};
struct c2 {};
struct c3 {};
struct c4 {};

struct beta {};
struct dirichlet {};
struct gendirichlet {};
struct gamma {};
struct slm {};
struct glm {};

struct title {};
struct selected {};
struct incpar {};
struct component {};
struct interval {};
struct io {};
struct param {};

struct pegtlstr {};

} // ctr::
} // quinoa::

#endif // QuinoaInputDeckTags_h
