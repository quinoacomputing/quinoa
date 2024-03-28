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
#include "TaggedTuple.hpp"

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
struct verbose { static std::string name() { return "verbose"; } };
struct nonblocking { static std::string name() { return "nonblocking"; } };
struct benchmark { static std::string name() { return "benchmark"; } };
struct lboff {};
struct feedback { static std::string name() { return "feedback"; } };
struct reorder { static std::string name() { return "reorder"; } };
struct lbfreq { static std::string name() { return "lbfreq"; } };
struct rsfreq { static std::string name() { return "rsfreq"; } };
struct pdf { static std::string name() { return "pdf"; } };
struct nchare {};
struct bounds {};
struct ordinary {};
struct central {};
struct binsize { static std::string name() { return "binsize"; } };
struct extent { static std::string name() { return "extent"; } };
struct tolderef { static std::string name() { return "tolderef"; } };
struct partitioner { static std::string name() { return "partitioner"; } };
struct partitioned { static std::string name() { return "partitioned"; } };
struct initpolicy { static std::string name() { return "initpolicy"; } };
struct coeffpolicy { static std::string name() { return "coeffpolicy"; } };
struct montecarlo {};
struct matched {};
struct compatibility {};
struct bndint {};
struct part { static std::string name() { return "part"; } };
struct particles { static std::string name() { return "particles"; } };
struct centroid {};
struct tty { static std::string name() { return "tty"; } };
struct dump {};
struct plot {};
struct glob {};
struct control { static std::string name() { return "control"; } };
struct stat { static std::string name() { return "stat"; } };
struct field { static std::string name() { return "field"; } };
struct surface { static std::string name() { return "surface"; } };
struct b { static std::string name() { return "b"; } };
struct S { static std::string name() { return "S"; } };
struct r { static std::string name() { return "r"; } };
struct c { static std::string name() { return "c"; } };
struct c0 { static std::string name() { return "c0"; } };
struct c1 {};
struct c2 {};
struct c3 { static std::string name() { return "c3"; } };
struct c4 { static std::string name() { return "c4"; } };
struct mean { static std::string name() { return "mean"; } };
struct cov { static std::string name() { return "cov"; } };
struct timescale { static std::string name() { return "timescale"; } };
struct virtualization {static std::string name() { return "virtualization"; }};
struct hydro {};
struct mix {};
struct frequency {};
struct mixrate {};
struct selected { static std::string name() { return "selected"; } };
struct discr { static std::string name() { return "discr"; } };
struct component { static std::string name() { return "component"; } };
struct iter { static std::string name() { return "iter"; } };
struct time { static std::string name() { return "time"; } };
struct range { static std::string name() { return "range"; } };
struct param { static std::string name() { return "param"; } };
struct init { static std::string name() { return "init"; } };
struct solve { static std::string name() { return "solve"; } };
struct chare { static std::string name() { return "chare"; } };
struct battery { static std::string name() { return "battery"; } };
struct generator {};
struct help { static std::string name() { return "help"; } };
struct helpctr { static std::string name() { return "helpctr"; } };
struct helpkw {};
struct cmdinfo {};
struct group { static std::string name() { return "group"; } };
struct esup {};
struct psup {};
struct gid {};
struct position_id { static std::string name() { return "position_id"; } };
struct velocity_id { static std::string name() { return "velocity_id"; } };
struct dissipation_id {
  static std::string name() { return "dissipation_id"; } };
struct mixmassfracbeta_id {
  static std::string name() { return "mixmassfracbeta_id"; } };
struct edge { static std::string name() { return "edge"; } };
struct com {};
struct queried {};
struct responded {};
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
struct update {};
struct ch {};
struct pe {};
struct it {};
struct node {};

DEFTAG(cmd);
DEFTAG(ctrinfo);

DEFTAG(title);
DEFTAG(nstep);
DEFTAG(term);
DEFTAG(t0);
DEFTAG(dt);
DEFTAG(cfl);
DEFTAG(ttyi);
DEFTAG(steady_state);
DEFTAG(residual);
DEFTAG(rescomp);
DEFTAG(partitioning);
DEFTAG(pelocal_reorder);
DEFTAG(operator_reorder);
DEFTAG(scheme);
DEFTAG(ndof);
DEFTAG(rdof);
DEFTAG(flux);
DEFTAG(limiter);
DEFTAG(cweight);
DEFTAG(shock_detector_coeff);
DEFTAG(accuracy_test);
DEFTAG(limsol_projection);
DEFTAG(fct);
DEFTAG(fctclip);
DEFTAG(fcteps);
DEFTAG(ctau);
DEFTAG(sysfct);
DEFTAG(sysfctvar);

DEFTAG(ncomp);
DEFTAG(pde);
DEFTAG(problem);
DEFTAG(transport);
DEFTAG(compflow);
DEFTAG(multimat);
DEFTAG(diffusivity);
DEFTAG(lambda);
DEFTAG(u0);
DEFTAG(alpha);
DEFTAG(beta);
DEFTAG(betax);
DEFTAG(betay);
DEFTAG(betaz);
DEFTAG(r0);
DEFTAG(p0);
DEFTAG(ce);
DEFTAG(kappa);
DEFTAG(nmat);
DEFTAG(prelax);
DEFTAG(prelax_timescale);
DEFTAG(intsharp);
DEFTAG(intsharp_param);
DEFTAG(depvar);
DEFTAG(sys);
DEFTAG(physics);

DEFTAG(material);
DEFTAG(eos);
DEFTAG(id);
DEFTAG(gamma);
DEFTAG(pstiff);
DEFTAG(w_gru);
DEFTAG(A_jwl);
DEFTAG(B_jwl);
DEFTAG(C_jwl);
DEFTAG(R1_jwl);
DEFTAG(R2_jwl);
DEFTAG(rho0_jwl);
DEFTAG(de_jwl);
DEFTAG(rhor_jwl);
DEFTAG(Tr_jwl);
DEFTAG(Pr_jwl);
DEFTAG(mu);
DEFTAG(cv);
DEFTAG(k);
DEFTAG(matidxmap);
DEFTAG(eosidx);
DEFTAG(matidx);
DEFTAG(solidx);

DEFTAG(field_output);
DEFTAG(iter_interval);
DEFTAG(time_interval);
DEFTAG(time_range);
DEFTAG(refined);
DEFTAG(filetype);
DEFTAG(sideset);
DEFTAG(outvar);
DEFTAG(diagnostics);
DEFTAG(error);
DEFTAG(format);
DEFTAG(precision);
DEFTAG(history_output);
DEFTAG(point);
DEFTAG(coord);

DEFTAG(ale);
DEFTAG(smoother);
DEFTAG(mesh_velocity);
DEFTAG(mesh_motion);
DEFTAG(meshforce);
DEFTAG(dvcfl);
DEFTAG(vortmult);
DEFTAG(maxit);
DEFTAG(tolerance);
DEFTAG(move);
DEFTAG(fntype);
DEFTAG(fn);

DEFTAG(amr);
DEFTAG(t0ref);
DEFTAG(dtref);
DEFTAG(dtref_uniform);
DEFTAG(dtfreq);
DEFTAG(maxlevels);
DEFTAG(initial);
DEFTAG(edgelist);
DEFTAG(coords);
DEFTAG(xminus);
DEFTAG(xplus);
DEFTAG(yminus);
DEFTAG(yplus);
DEFTAG(zminus);
DEFTAG(zplus);
DEFTAG(refvar);
DEFTAG(tol_refine);
DEFTAG(tol_derefine);
DEFTAG(pref);
DEFTAG(indicator);
DEFTAG(ndofmax);
DEFTAG(tolref);

DEFTAG(bc);
DEFTAG(mesh);
DEFTAG(dirichlet);
DEFTAG(symmetry);
DEFTAG(inlet);
DEFTAG(outlet);
DEFTAG(farfield);
DEFTAG(extrapolate);
DEFTAG(stag_point);
DEFTAG(sponge);
DEFTAG(vparam);
DEFTAG(pparam);
DEFTAG(radius);
DEFTAG(timedep);

DEFTAG(ic);
DEFTAG(materialid);
DEFTAG(pressure);
DEFTAG(temperature);
DEFTAG(velocity);
DEFTAG(box);
DEFTAG(meshblock);
DEFTAG(blockid);
DEFTAG(volume);
DEFTAG(mass);
DEFTAG(density);
DEFTAG(energy);
DEFTAG(energy_content);
DEFTAG(xmin);
DEFTAG(xmax);
DEFTAG(ymin);
DEFTAG(ymax);
DEFTAG(zmin);
DEFTAG(zmax);
DEFTAG(orientation);
DEFTAG(initiate);
DEFTAG(init_time);
DEFTAG(front_width);
DEFTAG(front_speed);

DEFTAG(filename);
DEFTAG(location);
DEFTAG(transfer);

struct BirthdaySpacings {};

} // tag::

#endif // Tags_h
