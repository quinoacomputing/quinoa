//******************************************************************************
/*!
  \file      src/Control/Quinoa/Tags.h
  \author    J. Bakosi
  \date      Mon 08 Sep 2014 07:45:31 PM MDT
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     Quinoa's tags
  \details   Quinoa's tags
*/
//******************************************************************************
#ifndef QuinoaTags_h
#define QuinoaTags_h

namespace quinoa {
namespace tag {

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
struct diffeq {};
struct virtualization {};

struct dirichlet {};
struct gendir {};
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
struct discr {};
struct component {};
struct interval {};
struct io {};
struct cmd {};
struct param {};
struct binsize {};

} // tag::
} // quinoa::

#endif // QuinoaTags_h
