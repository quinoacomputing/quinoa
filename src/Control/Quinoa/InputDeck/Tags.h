//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Tags.h
  \author    J. Bakosi
  \date      Thu Oct  3 10:15:49 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck tags
  \details   Quinoa's input dect tags
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
struct output {};
struct stat {};

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
