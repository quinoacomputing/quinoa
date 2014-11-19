//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Keywords.h
  \author    J. Bakosi
  \date      Tue 28 Oct 2014 09:25:33 PM MDT
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     Quinoa's input deck keywords
  \details   All keywords recognized by Quinoa's input deck parser. The keywords
  are defined by specializing struct 'keyword', defined in Control/Keyword.h.
  Introducing a new keyword requires a more human readable (but still short)
  name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef QuinoaInputDeckKeywords_h
#define QuinoaInputDeckKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes, such as *Keywords.h, below (if any) to make sure they
//! get included in the correct namespace and not polluting the global one.
#define Keywords

#include <Keyword.h>
#include <SharedKeywords.h>

namespace quinoa {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using tk::kw::keyword;
using tk::kw::undefined_info;

// Keyword 'hommix'
struct hommix_info {
  static const char* name() { return "Homogeneous material mixing"; }
  static const char* help() { return
    "Physics option, 'hommix', is short for homogeneous material mixing. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test material mixing models independent, i.e., decoupled from other "
    "equations. Only a set of scalar equations are advanced which can be "
    "coupled to each other. The keyword 'hommix' introduces the hommix ... end "
    "block, selecting and describing the parameters of the mixing model(s). "
    "The common physics keywords are recognized.";
  }
};
using hommix = keyword<hommix_info, h,o,m,m,i,x>;

// Keyword 'homhydro'
struct homhydro_info {
  static const char* name() { return "Homogeneous hydrodynamics"; }
  static const char* help() { return
    "Physics option, 'homhydro', is short for homogeneous hydrodynamics. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test hydrodynamics models independent of, i.e., decoupled from other "
    "equations. Only a set of momentum equations are advanced whose components "
    "can be coupled to each other. The keyword 'homhydro' introduces the "
    "homhydro ... end block, selecting and describing the parameters of the "
    "hydrodynamics model(s). The common physics keywords are recognized.";
  }
};
using homhydro = keyword<homhydro_info,  h,o,m,h,y,d,r,o >;

// Keyword 'homrt'
struct homrt_info {
  static const char* name() { return "Homogeneous Rayleigh-Taylor"; }
  static const char* help() { return
    "Physics option, 'homrt', is short for homogeneous Rayleigh-Taylor. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test hydrodynamics models for variable-density hydrodynamics and coupled "
    "material mixing, independent, i.e., decoupled from other equations. Only "
    "a set of mass and momentum conservation equations are advanced whose "
    "components can be coupled to each other. The keyword 'homrt' introduces "
    "the homrt ... end block, selecting and describing the parameters of the "
    "mass and hydrodynamics model(s). The common physics keywords are "
    "recognized.";
  }
};
using homrt = keyword<homrt_info,  h,o,m,r,t >;

// Keyword 'spinsflow'
struct spinsflow_info {
  static const char* name() { return
    "Standalone-particle incompressible Navier-Stokes flow";
  }
  static const char* help() { return
    "Physics option, 'spinsflow', is short for standalone-particle "
    "incompressible Navier-Stokes flow. It is a physics option intended for "
    "inhomogeneous constant-density flow. The transport equations solved are "
    "the momentum and optionally, energy, and a set of scalars. The "
    "divergence-constraint is enforced by solving a Poisson equation and "
    "projection scheme. The keyword 'spinsflow' introduces the spinsflow ... "
    "end block, selecting and describing the parameters of the above transport "
    "equations and their models. The common physics keywords are recognized.";
  }
};
using spinsflow = keyword<spinsflow_info,  s,p,i,n,s,f,l,o,w >;

// Keyword 'testsde'
struct testsde_info {
  static const char* name() {
    return "Differential equations testbed"; }
  static const char* help() { return
    "Test ordinary and stochastic differential equations.";
  }
};
using testsde = keyword<testsde_info,  t,e,s,t,s,d,e >;

// Keyword 'dirichlet'
struct dirichlet_info {
  static const char* name() {
    return "Dirichlet SDE"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "the Dirichlet distribution. For more details, see "
    "http://dx.doi.org/10.1155/2013/842981.";
  }
};
using dirichlet = keyword<dirichlet_info,  d,i,r,i,c,h,l,e,t >;

// Keyword 'generalized_dirichlet'
struct gendir_info {
  static const char* name() {
    return "Generalized Dirichlet SDE"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "Lochner's generalized Dirichlet distribution. For more details, see "
    "http://dx.doi.org/10.1063/1.4822416.";
  }
};
using gendir =
  keyword< gendir_info, g,e,n,e,r,a,l,i,z,e,d,'_',d,i,r,i,c,h,l,e,t >;

// Keyword 'wright_fisher'
struct wrightfisher_info {
  static const char* name() {
    return "Wright-Fisher SDE"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "the Dirichlet distribution. For more details, see "
    "http://www.sciencedirect.com/science/article/pii/S0040580912001013.";
  }
};
using wrightfisher =
  keyword< wrightfisher_info, w,r,i,g,h,t,'_',f,i,s,h,e,r >;

// Keyword 'skewnormal'
struct skewnormal_info {
  static const char* name() {
    return "Skew-normal SDE"; }
  static const char* help() { return
    "A single-variate stochastic differential equation whose invariant is the "
    "skew-normal distribution.";
  }
};
using skewnormal = keyword< skewnormal_info,  s,k,e,w,'-',n,o,r,m,a,l >;

// Keyword 'beta'
struct beta_info {
  static const char* name() {
    return "Beta SDE"; }
  static const char* help() { return
    "A single-variate stochastic differential equation whose invariant is the "
    "beta distribution.";
  }
};
using beta = keyword< beta_info,  b,e,t,a >;

// Keyword 'gamma'
struct gamma_info {
  static const char* name() {
    return "Gamma SDE"; }
  static const char* help() { return
    "A single-variate stochastic differential equation whose invariant is the "
    "gamma distribution.";
  }
};
using gamma = keyword< gamma_info,  g,a,m,m,a >;

// Keyword 'ornstein_uhlenbeck'
struct ornstein_uhlenbeck_info {
  static const char* name() {
    return "Ornstein-Uhlenbeck SDE"; }
  static const char* help() { return
    "A single-variate stochastic differential equation whose invariant is the "
    "normal distribution.";
  }
};
using ornstein_uhlenbeck =
  keyword< ornstein_uhlenbeck_info,  o,r,n,s,t,e,i,n,'-',u,h,l,e,n,b,e,c,k >;

// Keyword 'lognormal'
struct lognormal_info {
  static const char* name() {
    return "Log-normal SDE"; }
  static const char* help() { return
    "A single-variate stochastic differential equation whose invariant is the "
    "log-normal distribution.";
  }
};
using lognormal = keyword< lognormal_info,  l,o,g,n,o,r,m,a,l >;

// Select position model:
//   * Insviscid model
using pos_inviscid = keyword<undefined_info,  p,o,s,'_',i,n,v,i,s,c,i,d >;
//   * Viscous model
using pos_viscous = keyword<undefined_info,  p,o,s,'_',v,i,s,c,o,u,s >;

// Select mass model:
//   * Beta model
using mass_beta = keyword<undefined_info,  m,a,s,s,'_',b,e,t,a >;

// Select hydrodynamics model:
//   * Simplified Langevin model
using hydro_slm = keyword<undefined_info,  h,y,d,r,o,'_',s,l,m >;
//   * Generalized Langevin model
using hydro_glm = keyword<undefined_info,  h,y,d,r,o,'_',g,l,m >;

// Keyword 'mix_iem'
struct mix_iem_info {
  static const char* name() { return "Interaction by exchange with the mean"; }
  static const char* help() { return
    "Material mix model, 'mix_iem', is short for interaction by exchange with "
    "the mean (IEM). It is a relaxation-type material mix model intended "
    "shear-driven flows.";
  }
};
using mix_iem = keyword<mix_iem_info,  m,i,x,'_',i,e,m >;

// Keyword 'mix_iecm'
struct mix_iecm_info {
  static const char* name() { return
    "Interaction by exchange with the conditional mean";
  }
  static const char* help() { return
    "Material mix model, 'mix_iecm', is short for interaction by exchange with "
    "the conditional mean (IECM). It is a relaxation-type material mix model "
    "intended shear-driven flows.";
  }
};
using mix_iecm = keyword<mix_iecm_info,  m,i,x,'_',i,e,c,m >;

// Keyword 'mix_dir'
struct mix_dir_info {
  static const char* name() { return "Dirichlet"; }
  static const char* help() { return
    "Material mix model, 'mix_dir', is short for Dirichlet. It is a material "
    "mix model that explicitly satisfies the unit-sum requirement for all "
    "statistical samples.";
  }
};
using mix_dir = keyword<mix_dir_info,  m,i,x,'_',d,i,r >;

// Keyword 'mix_gendir'
struct mix_gendir_info {
  static const char* name() { return "Generalized Dirichlet"; }
  static const char* help() { return
    "Material mix model, 'mix_gendir', is short for Lochner's generalized "
    "Dirichlet. It is a material mix model that explicitly satisfies the "
    "unit-sum requirement for all statistical samples.";
  }
};
using mix_gendir = keyword<mix_gendir_info,  m,i,x,'_',g,e,n,d,i,r >;

// Select material mix rate model:
//   * Gamma distribution model
using mixrate_gamma = keyword<undefined_info,  m,i,x,r,a,t,e,'_',g,a,m,m,a >;

// Select turbulence frequency model:
//   * Gamma distribution model
using freq_gamma = keyword<undefined_info,  f,r,e,q,'_',g,a,m,m,a >;

// Number of time steps to take
using nstep = keyword<undefined_info,  n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword<undefined_info,  t, e, r, m >;

// Size of time step
using dt = keyword<undefined_info,  d,t >;

// Start of position model specification block
using position = keyword<undefined_info,  p,o,s,i,t,i,o,n >;

// Start of hydrodynamics model specification block
using hydro = keyword<undefined_info,  h,y,d,r,o >;

// Start of material mix model specification block
using mix = keyword<undefined_info,  m,i,x >;

// Number of particle position components
using nposition = keyword<undefined_info,  n,p,o,s,i,t,i,o,n >;
// Number of particle density components
using ndensity = keyword<undefined_info,  n,d,e,n,s,i,t,y >;
// Number of particle velocity components
using nvelocity = keyword<undefined_info,  n,v,e,l,o,c,i,t,y >;
// Number of particle scalar components
using nscalar = keyword<undefined_info,  n,s,c,a,l,a,r >;
// Number of particle turbulence frequency components
using nfreq = keyword<undefined_info,  n,f,r,e,q >;

// Number of components
using ncomp = keyword<undefined_info,  n,c,o,m,p >;

// Wright-Fisher SDE parameters
using sde_omega = keyword<undefined_info, o,m,e,g,a >;

// Dirichlet and generalized Dirichlet parameters
using sde_b = keyword<undefined_info,  b >;
using sde_S = keyword<undefined_info,  S >;
using sde_kappa = keyword<undefined_info,  k,a,p,p,a >;
using sde_c = keyword<undefined_info,  c >;

// Ornstein-Uhlenbeck SDE parameters
using sde_sigma = keyword<undefined_info, s,i,g,m,a >;
using sde_theta = keyword<undefined_info, t,h,e,t,a >;
using sde_mu = keyword<undefined_info, m,u >;

// Langevin model parameters
using SLM_C0 = keyword<undefined_info,  C,'0' >;

// Gamma frequency model parameters
using freq_gamma_C1 = keyword<undefined_info,  C,'1' >;
using freq_gamma_C2 = keyword<undefined_info,  C,'2' >;
using freq_gamma_C3 = keyword<undefined_info,  C,'3' >;
using freq_gamma_C4 = keyword<undefined_info,  C,'4' >;

// Beta model parameters
using Beta_At = keyword<undefined_info,  A,t >;

// Quantities
using transported_scalar = keyword<undefined_info,  Y >;
using transported_scalar_fluctuation = keyword<undefined_info,  y >;

using velocity_x = keyword<undefined_info,  U >;
using velocity_fluctuation_x = keyword<undefined_info,  u >;
using velocity_y = keyword<undefined_info,  V >;
using velocity_fluctuation_y = keyword<undefined_info,  v >;
using velocity_z = keyword<undefined_info,  W >;
using velocity_fluctuation_z = keyword<undefined_info,  w >;

using pressure = keyword<undefined_info,  P >;
using pressure_fluctuation = keyword<undefined_info,  p >;
  
using density = keyword<undefined_info,  R >;
using density_fluctuation = keyword<undefined_info,  r >;
  
// Total number of particles
using npar = keyword<undefined_info,  n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword<undefined_info,  t,t,y,i >;

// Dump (restart file) output interval
using dmpi = keyword<undefined_info,  d,m,p,i >;

// Glob output interval
using glbi = keyword<undefined_info,  g,l,b,i >;

// Output interval
using interval = keyword<undefined_info,  i,n,t,e,r,v,a,l >;

// Statistics
using statistics = keyword<undefined_info,  s,t,a,t,i,s,t,i,c,s >;

// PDFs
using pdfs = keyword<undefined_info, p,d,f,s >;

// Keyword 'filetype': selected PDF output file type
struct filetype_info {
  static const char* name() { return "PDF output file type"; }
  static const char* help() { return "PDF output file type."; }
};
using pdf_filetype = keyword< filetype_info, f,i,l,e,t,y,p,e >;

// Keyword 'txt': PDF output file type: txt
struct txt_info {
  static const char* name() { return "text"; }
  static const char* help() { return "PDF text output file."; }
};
using txt = keyword< txt_info, t,x,t >;

// Keyword 'gmshtxt': PDF output file type: text gmsh
struct gmshtxt_info {
  static const char* name() { return "gmshtxt"; }
  static const char* help() { return "PDF GMSH txt file output."; }
};
using gmshtxt = keyword< gmshtxt_info, g,m,s,h,t,x,t >;

// Keyword 'gmshbin': PDF output file type: binary gmsh
struct gmshbin_info {
  static const char* name() { return "gmshbin"; }
  static const char* help() { return "PDF GMSH binary file output."; }
};
using gmshbin = keyword< gmshbin_info, g,m,s,h,b,i,n >;

// Keyword 'exodusii': PDF output file type: Exodus II
struct exodusii_info {
  static const char* name() { return "exodusii"; }
  static const char* help() { return "PDF Exodus II binary file output."; }
};
using exodusii = keyword< exodusii_info, e,x,o,d,u,s,i,i >;

// Keyword 'policy': selected PDF output file policy
struct policy_info {
  static const char* name() { return "PDF output file policy"; }
  static const char* help() { return "PDF output file policy."; }
};
using pdf_policy = keyword< policy_info, p,o,l,i,c,y >;

// Keyword 'overwrite': PDF output file policy: overwrite
struct overwrite_info {
  static const char* name() { return "overwrite"; }
  static const char* help() { return "PDF output overwrites the same file "
  "containing a single time step."; }
};
using overwrite = keyword< overwrite_info, o,v,e,r,w,r,i,t,e >;

// Keyword 'multiple': PDF output file policy: multiple
struct multiple_info {
  static const char* name() { return "multiple"; }
  static const char* help() { return "PDF output creates new file for each "
  "time step."; }
};
using multiple = keyword< multiple_info, m,u,l,t,i,p,l,e >;

// Keyword 'evolution': PDF output file policy: evolution
struct evolution_info {
  static const char* name() { return "evolution"; }
  static const char* help() { return "PDF output appends new time step to the "
  "same file for each yielding a time evolution of data."; }
};
using evolution = keyword< evolution_info, e,v,o,l,u,t,i,o,n >;

// Keyword 'txt_float_format': PDF txt output floating-point format
struct txt_float_format_info {
  static const char* name() { return "format"; }
  static const char* help() { return "Formatting for floating-point output."; }
};
using txt_float_format = keyword< txt_float_format_info, f,o,r,m,a,t >;

// Keyword 'txt_float_default': PDF txt output floating-point format: default
struct txt_float_default_info {
  static const char* name() { return "default"; }
  static const char* help() { return "Use default formatting for "
  "floating-point output."; }
};
using txt_float_default = keyword< txt_float_default_info, d,e,f,a,u,l,t >;

// Keyword 'txt_float_scientific': PDF txt output floating-point format:
// scientific
struct txt_float_scientific_info {
  static const char* name() { return "scientific"; }
  static const char* help() { return "Use scientific formatting for "
  "floating-point output."; }
};
using txt_float_scientific =
  keyword< txt_float_scientific_info, s,c,i,e,n,t,i,f,i,c >;

// Keyword 'fixed': PDF txt output floating-point format: fixed
struct txt_float_fixed_info {
  static const char* name() { return "fixed"; }
  static const char* help() { return "Use fixed formatting for "
  "floating-point output."; }
};
using txt_float_fixed = keyword< txt_float_fixed_info, f,i,x,e,d >;

// Keyword 'precision': precision in digits for text output of numbers
struct precision_info {
  static const char* name() { return "precision"; }
  static const char* help() { return "Precision in digits for text output.."; }
};
using precision = keyword< precision_info, p,r,e,c,i,s,i,o,n >;

// Keyword 'centering': selected PDF output file centering
struct centering_info {
  static const char* name() { return "PDF output file centering"; }
  static const char* help() { return "PDF output file centering."; }
};
using pdf_centering = keyword< centering_info, c,e,n,t,e,r,i,n,g >;

// Keyword 'elem': PDF output file centering: elem
struct elem_info {
  static const char* name() { return "elem"; }
  static const char* help() { return "PDF output bins at elem centers."; }
};
using elem = keyword< elem_info, e,l,e,m >;

// Keyword 'node': PDF output file centering: node
struct node_info {
  static const char* name() { return "node"; }
  static const char* help() { return "PDF output bins at nodes."; }
};
using node = keyword< node_info, n,o,d,e >;

// RNG block
using rngs = keyword<undefined_info,  r,n,g,s >;

// RNG
using rng = keyword<undefined_info,  r,n,g >;

// Keyword 'init'
using init = keyword<undefined_info,  i,n,i,t >;

// Keyword 'coeff'
using coeff = keyword<undefined_info,  c,o,e,f,f >;

// Keyword 'raw': raw initialization policy
struct raw_info {
  static const char* name() { return "Raw"; }
  static const char* help() { return "Raw initialization policy."; }
};
using raw = keyword<raw_info, r,a,w >;

// Keyword 'raw': zero initialization policy
struct zero_info {
  static const char* name() { return "Zero"; }
  static const char* help() { return "Zero initialization policy."; }
};
using zero = keyword<zero_info, z,e,r,o >;

// Keyword 'constant': constant coefficients policy
struct constant_info {
  static const char* name() { return "Constant"; }
  static const char* help() { return "Constant coefficients policy."; }
};
using constant = keyword<constant_info, c,o,n,s,t,a,n,t >;

// Keyword 'depvar': dependent variable in equations
struct depvar_info {
  static const char* name() { return "Dependent variable"; }
  static const char* help() { return "Dependent variable in differential "
                                     "equations."; }
};
using depvar = keyword< depvar_info, d,e,p,v,a,r >;


} // kw::
} // quinoa::

#undef Keywords

#endif // QuinoaInputDeckKeywords_h
