//******************************************************************************
/*!
  \file      src/Control/QuinoaControlTags.h
  \author    J. Bakosi
  \date      Thu Sep 12 10:35:46 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's control tags
  \details   All control tags used to build a nested tagged tagged tuple
*/
//******************************************************************************
#ifndef QuinoaControlTags_h
#define QuinoaControlTags_h

namespace quinoa {

namespace control {

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

struct ctr {};
struct input {};
struct output {};
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

} // namespace control

} // namespace quinoa

#endif // QuinoaControlTags_h
