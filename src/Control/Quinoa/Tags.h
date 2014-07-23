//******************************************************************************
/*!
  \file      src/Control/Quinoa/Tags.h
  \author    J. Bakosi
  \date      Wed 23 Jul 2014 07:16:24 AM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's tags
  \details   Quinoa's tags
*/
//******************************************************************************
#ifndef QuinoaTags_h
#define QuinoaTags_h

namespace quinoa {
namespace tag {

struct geometry {};
struct montecarlo {};
struct nstep {};
struct term {};
struct dt {};
struct npar {};
struct ncomp {};
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

struct dirichlet { enum { id = 0 }; };
struct gendir { enum { id = 1 }; };
struct beta {};
struct gamma {};
struct slm {};
struct glm {};
struct ou {};
struct lognormal {};
struct skewnormal {};
struct position {};
struct mass {};
struct hydro {};
struct energy {};
struct mix {};
struct frequency {};
struct mixrate {};

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
