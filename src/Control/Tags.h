//******************************************************************************
/*!
  \file      src/Control/Tags.h
  \author    J. Bakosi
  \date      Fri 17 Apr 2015 08:53:47 AM MDT
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Tags
  \details   Tags are unique types, used for metaprogramming.
*/
//******************************************************************************
#ifndef Tags_h
#define Tags_h

//! Tags used as unique-type labels for compile-time code-generation
namespace tag {

struct io {};
struct input {};
struct output {};
struct seed {};
struct uniform_method {};
struct gaussian_method {};
struct rng {};
struct rngmkl {};
struct rngsse {};
struct seqlen {};
struct verbose {};
struct error {};
struct pdf {};
struct ordpdf {};
struct cenpdf {};
struct pdffiletype {};
struct pdfpolicy {};
struct pdfctr {};
struct pdfnames {};
struct flformat {};
struct prec {};
struct ordinary {};
struct central {};
struct binsize {};
struct extent {};
struct dirichlet {};
struct gendir {};
struct wrightfisher {};
struct beta {};
struct numfracbeta {};
struct massfracbeta {};
struct mixnumfracbeta {};
struct mixmassfracbeta {};
struct gamma {};
struct spike {};
struct diffeq {};
struct initpolicy {};
struct coeffpolicy {};
struct montecarlo {};
struct nstep {};
struct term {};
struct dt {};
struct npar {};
struct ncomp {};
struct tty {};
struct dump {};
struct plot {};
struct glob {};
struct control {};
struct stat {};
struct atwood {};
struct b {};
struct S {};
struct kappa {};
struct bprime {};
struct kappaprime {};
struct rho2 {};
struct rcomma {};
struct r {};
struct c {};
struct c0 {};
struct c1 {};
struct c2 {};
struct c3 {};
struct c4 {};
struct lambda {};
struct sigmasq {};
struct theta {};
struct mu {};
struct timescale {};
struct depvar {};
struct virtualization {};
struct omega {};
struct slm {};
struct glm {};
struct diagou {};
struct ou {};
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
struct cmd {};
struct param {};
struct init {};
struct chare {};
struct battery {};
struct generator {};
struct help {};
struct helpctr {};
struct helpkw {};
struct cmdinfo {};
struct ctrinfo {};
struct group {};
struct esup {};
struct psup {};
struct owner {};

struct BirthdaySpacings {};
struct Collision {};
struct RandomWalk1 {};
struct Gap {};
struct SimplePoker {};
struct CouponCollector {};
struct MaxOft {};
struct WeightDistrib {};
struct MatrixRank {};
struct HammingIndep {};
struct SerialOver {};
struct CollisionOver {};
struct ClosePairs {};
struct ClosePairsBitMatch {};
struct Run {};
struct Permutation {};
struct CollisionPermut {};
struct SampleProd {};
struct SampleMean {};
struct SampleCorr {};
struct AppearanceSpacings {};
struct SumCollector {};
struct Savir2 {};
struct GCD {};
struct LinearComp {};
struct LempelZiv {};
struct Fourier3 {};
struct LongestHeadRun {};
struct PeriodsInStrings {};
struct HammingWeight2 {};
struct HammingCorr {};
struct StringRun {};
struct AutoCorr {};

} // tag::

#endif // Tags_h
