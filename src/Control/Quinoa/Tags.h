//******************************************************************************
/*!
  \file      src/Control/Quinoa/Tags.h
  \author    J. Bakosi
  \date      Fri 10 Oct 2014 02:50:31 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
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
struct theta {};
struct mu {};
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
struct extent {};
struct pdffiletype {};
struct pdfpolicy {};
struct pdfctr {};
struct pdfnames {};

struct float_format {};
struct precision {};
struct init {};
struct ordinary {};
struct central {};
struct chare {};

} // tag::
} // quinoa::

#endif // QuinoaTags_h
