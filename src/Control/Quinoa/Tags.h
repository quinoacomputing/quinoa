//******************************************************************************
/*!
  \file      src/Control/Quinoa/Tags.h
  \author    J. Bakosi
  \date      Wed 19 Feb 2014 05:30:15 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's input deck tags
  \details   Quinoa's input dect tags
*/
//******************************************************************************
#ifndef QuinoaTags_h
#define QuinoaTags_h

namespace quinoa {
namespace tag {

struct geometry {};
struct montecarlo {};
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
struct ncomp {};
struct ndirichlet {};
struct ngendir {};
struct nou {};
struct nlognormal {};
struct nskewnormal {};
struct ngamma {};
struct nbeta {};

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
struct lambda {};
struct sigma {};
struct timescale {};

struct initpolicy {};
struct coeffpolicy {};

struct depvar {};

struct sde {};
struct beta {};
struct dirichlet {};
struct gendir {};
struct gamma {};
struct slm {};
struct glm {};
struct ou {};
struct lognormal {};
struct skewnormal {};

struct title {};
struct selected {};
struct incpar {};
struct component {};
struct interval {};
struct io {};
struct cmd {};
struct param {};

} // tag::
} // quinoa::

#endif // QuinoaTags_h
