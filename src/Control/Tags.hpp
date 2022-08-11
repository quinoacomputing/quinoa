// *****************************************************************************
/*!
  \file      src/Control/Tags.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
  \brief     Tags
  \details   Tags are unique types, used for metaprogramming.
*/
// *****************************************************************************
#ifndef Tags_h
#define Tags_h

#include <string>

//! Tags used as unique-type labels for compile-time code-generation
namespace tag {

struct low {};
struct high {};
struct io { static std::string name() { return "io"; } };
struct quiescence { static std::string name() { return "quiescence"; } };
struct trace { static std::string name() { return "trace"; } };
struct version { static std::string name() { return "version"; } };
struct license { static std::string name() { return "license"; } };
struct input { static std::string name() { return "input"; } };
struct output { static std::string name() { return "output"; } };
struct screen { static std::string name() { return "screen"; } };
struct restart { static std::string name() { return "restart"; } };
struct nrestart { static std::string name() { return "nrestart"; } };
struct diag { static std::string name() { return "diag"; } };
struct history { static std::string name() { return "history"; } };
struct evalLB {};
struct seed { static std::string name() { return "seed"; } };
struct uniform_method {
  static std::string name() { return "uniform_method"; } };
struct gaussian_method {
  static std::string name() { return "gaussian_method"; } };
struct gaussianmv_method {
  static std::string name() { return "gaussianmv_method"; } };
struct gaussian { static std::string name() { return "gaussian"; } };
struct jointgaussian { static std::string name() { return "jointgaussian"; } };
struct beta_method { static std::string name() { return "beta_method"; } };
struct gamma_method { static std::string name() { return "gamma_method"; } };
struct rng { static std::string name() { return "rng"; } };
struct xminus { static std::string name() { return "xminus"; } };
struct xplus { static std::string name() { return "xplus"; } };
struct yminus { static std::string name() { return "yminus"; } };
struct yplus { static std::string name() { return "yplus"; } };
struct zminus { static std::string name() { return "zminus"; } };
struct zplus { static std::string name() { return "zplus"; } };
struct rngmkl { static std::string name() { return "rngmkl"; } };
struct rngsse { static std::string name() { return "rngsse"; } };
struct rng123 { static std::string name() { return "rng123"; } };
struct seqlen { static std::string name() { return "seqlen"; } };
struct verbose { static std::string name() { return "verbose"; } };
struct nonblocking { static std::string name() { return "nonblocking"; } };
struct benchmark { static std::string name() { return "benchmark"; } };
struct lboff {};
struct feedback { static std::string name() { return "feedback"; } };
struct reorder { static std::string name() { return "reorder"; } };
struct pelocal_reorder {
  static std::string name() { return "pelocal_reorder"; } };
struct operator_reorder {
  static std::string name() { return "operator_reorder"; } };
struct steady_state {
  static std::string name() { return "steady_state"; } };
struct residual { static std::string name() { return "residual"; } };
struct error { static std::string name() { return "error"; } };
struct lbfreq { static std::string name() { return "lbfreq"; } };
struct rsfreq { static std::string name() { return "rsfreq"; } };
struct dtfreq { static std::string name() { return "dtfreq"; } };
struct pdf { static std::string name() { return "pdf"; } };
struct ordpdf {};
struct cenpdf {};
struct nchare {};
struct bounds {};
struct meshvelocity { static std::string name() { return "meshvelocity"; } };
struct smoother { static std::string name() { return "smoother"; } };
struct meshforce { static std::string name() { return "meshforce"; } };
struct mesh_motion{ static std::string name() { return "mesh_motion"; } };
struct vortmult { static std::string name() { return "vortmult"; } };
struct maxit { static std::string name() { return "maxit"; } };
struct tolerance { static std::string name() { return "tolerace"; } };
struct mesh { static std::string name() { return "mesh"; } };
struct couple { static std::string name() { return "couple"; } };
struct transfer { static std::string name() { return "transfer"; } };
struct filetype { static std::string name() { return "filetype"; } };
struct filename { static std::string name() { return "filename"; } };
struct location { static std::string name() { return "location"; } };
struct orientation { static std::string name() { return "orientation"; } };
struct reference { static std::string name() { return "reference"; } };
struct pdfpolicy { static std::string name() { return "pdfpolicy"; } };
struct pdfctr { static std::string name() { return "pdfctr"; } };
struct pdfnames { static std::string name() { return "pdfnames"; } };
struct flformat { static std::string name() { return "flformat"; } };
struct prec { static std::string name() { return "precision"; } };
struct ordinary {};
struct central {};
struct binsize { static std::string name() { return "binsize"; } };
struct extent { static std::string name() { return "extent"; } };
struct dirichlet { static std::string name() { return "dirichlet"; } };
struct mixdirichlet { static std::string name() { return "mixdirichlet"; } };
struct gendir { static std::string name() { return "gendir"; } };
struct wrightfisher { static std::string name() { return "wrightfisher"; } };
struct p0 { static std::string name() { return "p0"; } };
struct beta { static std::string name() { return "beta"; } };
struct betax { static std::string name() { return "betax"; } };
struct betay { static std::string name() { return "betay"; } };
struct betaz { static std::string name() { return "betaz"; } };
struct r0 { static std::string name() { return "r0"; } };
struct ce { static std::string name() { return "ce"; } };
struct numfracbeta { static std::string name() { return "numfracbeta"; } };
struct massfracbeta { static std::string name() { return "massfracbeta"; } };
struct mixnumfracbeta {
  static std::string name() { return "mixnumfracbeta"; } };
struct mixmassfracbeta {
  static std::string name() { return "mixmassfracbeta"; } };
struct alpha { static std::string name() { return "alpha"; } };
struct gamma { static std::string name() { return "gamma"; } };
struct pstiff { static std::string name() { return "pstiff"; } };
struct spike { static std::string name() { return "spike"; } };
struct betapdf { static std::string name() { return "betapdf"; } };
struct hydrotimescales {
  static std::string name() { return "hydrotimescales"; } };
struct hydroproductions {
  static std::string name() { return "hydroproductions"; } };
struct diffeq { static std::string name() { return "diffeq"; } };
struct pde { static std::string name() { return "pde"; } };
struct pref { static std::string name() { return "pref"; } };
struct tolref { static std::string name() { return "tolref"; } };
struct ndofmax { static std::string name() { return "ndofmax"; } };
struct indicator{ static std::string name() { return "indicator"; } };
struct shock_indicator {static std::string name() {return "shock_indicator";}};
struct mesh_size { static std::string name() { return "mesh_size"; } };
struct coeff { static std::string name() { return "coeff"; } };
struct amr { static std::string name() { return "amr"; } };
struct ale { static std::string name() { return "ale"; } };
struct move { static std::string name() { return "move"; } };
struct tolderef { static std::string name() { return "tolderef"; } };
struct t0ref { static std::string name() { return "t0ref"; } };
struct dtref { static std::string name() { return "dtref"; } };
struct dtref_uniform { static std::string name() { return "dtref_uniform"; } };
struct maxlevels { static std::string name() { return "maxlevels"; } };
struct partitioner { static std::string name() { return "partitioner"; } };
struct partitioned { static std::string name() { return "partitioned"; } };
struct scheme { static std::string name() { return "scheme"; } };
struct initpolicy { static std::string name() { return "initpolicy"; } };
struct coeffpolicy { static std::string name() { return "coeffpolicy"; } };
struct montecarlo {};
struct nstep { static std::string name() { return "nstep"; } };
struct term { static std::string name() { return "term"; } };
struct t0 { static std::string name() { return "t0"; } };
struct dt { static std::string name() { return "dt"; } };
struct cfl { static std::string name() { return "cfl"; } };
struct dvcfl { static std::string name() { return "dvcfl"; } };
struct fct { static std::string name() { return "fct"; } };
struct fctclip { static std::string name() { return "fctclip"; } };
struct sys { static std::string name() { return "sys"; } };
struct sysfct { static std::string name() { return "sysfct"; } };
struct sysfctvar { static std::string name() { return "sysfctvar"; } };
struct fcteps { static std::string name() { return "fcteps"; } };
struct ctau { static std::string name() { return "ctau"; } };
struct npar { static std::string name() { return "npar"; } };
struct refined { static std::string name() { return "refined output"; } };
struct matched {};
struct compatibility {};
struct bndint {};
struct part { static std::string name() { return "part"; } };
struct particles { static std::string name() { return "particles"; } };
struct centroid {};
struct ncomp { static std::string name() { return "ncomp"; } };
struct nmat { static std::string name() { return "nmat"; } };
struct prelax { static std::string name() { return "prelax"; } };
struct prelax_timescale {
  static std::string name() { return "prelax_timescales"; } };
struct intsharp { static std::string name() { return "intsharp"; } };
struct intsharp_param {
  static std::string name() { return "intsharp_param"; } };
struct tty { static std::string name() { return "tty"; } };
struct dump {};
struct plot {};
struct glob {};
struct control { static std::string name() { return "control"; } };
struct stat { static std::string name() { return "stat"; } };
struct field { static std::string name() { return "field"; } };
struct surface { static std::string name() { return "surface"; } };
struct atwood {};
struct b { static std::string name() { return "b"; } };
struct S { static std::string name() { return "S"; } };
struct kappa { static std::string name() { return "kappa"; } };
struct bprime { static std::string name() { return "bprime"; } };
struct kappaprime { static std::string name() { return "kappaprime"; } };
struct rho2 { static std::string name() { return "rho2"; } };
struct rho { static std::string name() { return "rho"; } };
struct mean_gradient { static std::string name() { return "mean_gradient"; } };
struct rcomma { static std::string name() { return "rcomma"; } };
struct r { static std::string name() { return "r"; } };
struct c { static std::string name() { return "c"; } };
struct c0 { static std::string name() { return "c0"; } };
struct gravity { static std::string name() { return "gravity"; } };
struct c1 {};
struct c2 {};
struct c3 { static std::string name() { return "c3"; } };
struct c4 { static std::string name() { return "c4"; } };
struct com1 { static std::string name() { return "com1"; } };
struct com2 { static std::string name() { return "com2"; } };
struct lambda { static std::string name() { return "lambda"; } };
struct sigmasq { static std::string name() { return "sigmasq"; } };
struct theta { static std::string name() { return "theta"; } };
struct mu { static std::string name() { return "mu"; } };
struct mean { static std::string name() { return "mean"; } };
struct cov { static std::string name() { return "cov"; } };
struct timescale { static std::string name() { return "timescale"; } };
struct depvar { static std::string name() { return "depvar"; } };
struct refvar { static std::string name() { return "refvar"; } };
struct virtualization {static std::string name() { return "virtualization"; }};
struct omega { static std::string name() { return "omega"; } };
struct slm { static std::string name() { return "slm"; } };
struct glm { static std::string name() { return "glm"; } };
struct diagou { static std::string name() { return "diagou"; } };
struct ou { static std::string name() { return "ou"; } };
struct skewnormal { static std::string name() { return "skewnormal"; } };
struct position { static std::string name() { return "position"; } };
struct dissipation { static std::string name() { return "dissipation"; } };
struct variant { static std::string name() { return "variant"; } };
struct normalization { static std::string name() { return "normalization"; } };
struct materialid { static std::string name() { return "materialid"; } };
struct matidxmap { static std::string name() { return "matidxmap"; } };
struct matidx { static std::string name() { return "matidx"; } };
struct eosidx { static std::string name() { return "eosidx"; } };
struct mass { static std::string name() { return "mass"; } };
struct hydro {};
struct mix {};
struct frequency {};
struct mixrate {};
struct title { static std::string name() { return "title"; } };
struct selected { static std::string name() { return "selected"; } };
struct discr { static std::string name() { return "discr"; } };
struct bc { static std::string name() { return "bc"; } };
struct farfield_pressure {
  static std::string name() { return "farfield_pressure"; } };
struct farfield_density {
  static std::string name() { return "farfield_density"; } };
struct farfield_velocity {
  static std::string name() { return "farfield_velocity"; } };
struct bcfarfield { static std::string name() { return "bcfarfield"; } };
struct component { static std::string name() { return "component"; } };
struct rescomp { static std::string name() { return "residual component"; } };
struct iter { static std::string name() { return "iter"; } };
struct time { static std::string name() { return "time"; } };
struct range { static std::string name() { return "range"; } };
struct cmd { static std::string name() { return "cmd"; } };
struct param { static std::string name() { return "param"; } };
struct init { static std::string name() { return "init"; } };
struct initiate { static std::string name() { return "initiate"; } };
struct solve { static std::string name() { return "solve"; } };
struct chare { static std::string name() { return "chare"; } };
struct battery { static std::string name() { return "battery"; } };
struct generator {};
struct help { static std::string name() { return "help"; } };
struct helpctr { static std::string name() { return "helpctr"; } };
struct helpkw {};
struct cmdinfo {};
struct ctrinfo {};
struct group { static std::string name() { return "group"; } };
struct esup {};
struct psup {};
struct gid {};
struct transport { static std::string name() { return "transport"; } };
struct compflow { static std::string name() { return "compflow"; } };
struct multimat { static std::string name() { return "multimat"; } };
struct problem { static std::string name() { return "problem"; } };
struct physics { static std::string name() { return "physics"; } };
struct diffusivity { static std::string name() { return "diffusivity"; } };
struct u0 { static std::string name() { return "u0"; } };
struct bcdir { static std::string name() { return "bcdir"; } };
struct bcsym { static std::string name() { return "bcsym"; } };
struct stag { static std::string name() { return "stag"; } };
struct skip { static std::string name() { return "skip"; } };
struct point { static std::string name() { return "point"; } };
struct radius { static std::string name() { return "radius"; } };
struct sideset { static std::string name() { return "sideset"; } };
struct sponge { static std::string name() { return "sponge"; } };
struct bcinlet { static std::string name() { return "bcinlet"; } };
struct bcoutlet { static std::string name() { return "bcoutlet"; } };
struct bcextrapolate { static std::string name() { return "bcextrapolate"; } };
struct bctimedep { static std::string name() { return "bctimedep"; } };
struct material { static std::string name() { return "material"; } };
struct eos { static std::string name() { return "eos"; } };
struct id { static std::string name() { return "id"; } };
struct position_id { static std::string name() { return "position_id"; } };
struct velocity_id { static std::string name() { return "velocity_id"; } };
struct dissipation_id {
  static std::string name() { return "dissipation_id"; } };
struct mixmassfracbeta_id {
  static std::string name() { return "mixmassfracbeta_id"; } };
struct edge { static std::string name() { return "edge"; } };
struct cv { static std::string name() { return "cv"; } };
struct k { static std::string name() { return "k"; } };
struct com {};
struct queried {};
struct responded {};
struct coord {};
struct refinserted {};
struct discinserted {};
struct disccreated {};
struct workinserted {};
struct distributed {};
struct load {};
struct bcast {};
struct elem {};
struct avecost {};
struct stdcost {};
struct flux { static std::string name() { return "flux"; } };
struct ndof{ static std::string name() { return "ndof"; } };
struct rdof{ static std::string name() { return "rdof"; } };
struct limiter { static std::string name() { return "limiter"; } };
struct accuracy_test { static std::string name() { return "accuracy_test"; } };
struct limsol_projection { static std::string name() { return "limsol_projection"; } };
struct shock_detection { static std::string name() { return "shock_detection"; } };
struct cweight { static std::string name() { return "cweight"; } };
struct update {};
struct ch {};
struct pe {};
struct it {};
struct fn { static std::string name() { return "fn"; } };
struct fntype { static std::string name() { return "fntype"; } };
struct node {};
struct ic { static std::string name() { return "ic"; } };
struct velocity { static std::string name() { return "velocity"; } };
struct density { static std::string name() { return "density"; } };
struct pressure { static std::string name() { return "pressure"; } };
struct energy { static std::string name() { return "energy"; } };
struct energy_content { static std::string name() { return "energy_content"; } };
struct temperature { static std::string name() { return "temperature"; } };
struct outvar { static std::string name() { return "outvar"; } };
struct box { static std::string name() { return "box"; } };
struct xmin { static std::string name() { return "xmin"; } };
struct xmax { static std::string name() { return "xmax"; } };
struct ymin { static std::string name() { return "ymin"; } };
struct ymax { static std::string name() { return "ymax"; } };
struct zmin { static std::string name() { return "zmin"; } };
struct zmax { static std::string name() { return "zmax"; } };
struct lua { static std::string name() { return "lua"; } };

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
