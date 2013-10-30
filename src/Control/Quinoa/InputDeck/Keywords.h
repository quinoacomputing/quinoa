//******************************************************************************
/*!
  \file      src/Control/Quinoa/InputDeck/Keywords.h
  \author    J. Bakosi
  \date      Wed Oct 30 06:45:14 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
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

namespace quinoa {
//! List of keywords the parser understands
namespace kw {

using namespace pegtl::ascii;
using namespace tk::kw;

// Include base keywords recognized by all input deck parsers
#include <BaseKeywords.h>

// Include Intel's MKL's RNG keywords
#include <MKLRNGKeywords.h>

// Keyword 'analytic_geometry'
struct analytic_geometry_info {
  static const char* name() { return "Analytic"; }
  static const char* help() { return
    "This option is used to define the beginning of an analytical geometry "
    "definition block. Analytical geometry definitions describe a 3D geometry "
    "using available primitives. Example:\n"
    "\tanalytical_geometry\n"
    "\t  box 0.0 0.0 0.0  1.1 1.2 1.3 end # opposite corners: x1 y1 z1 x2 y2 "
    "z2\n"
    "\tend";
  }
};
using analytic_geometry =
  keyword<analytic_geometry_info, a,n,a,l,y,t,i,c,'_',g,e,o,m,e,t,r,y>;

// Keyword 'discrete_geometry'
struct discrete_geometry_info {
  static const char* name() { return "Discrete"; }
  static const char* help() { return
    "This option is used to define the beginning of a discrete geometry "
    "definition block. Discrete geometry definitions are used to specify "
    "various parameters that influence the generation of point clouds within a "
    "closed geometry specified by discrete surfaces. Example:\n"
    "\tdiscrete_geometry\n"
    "\t  ...\n"
    "\tend";
  }
};
using discrete_geometry =
  keyword<discrete_geometry_info, d,i,s,c,r,e,t,e,'_',g,e,o,m,e,t,r,y>;

// Keyword 'brick'
struct brick_info {
  static const char* name() { return "Brick primitive"; }
  static const char* help() { return
    "A brick primitive is one of the available geometric primitives that can "
    "be used to define a 3D geometry in the analytical way, i.e., inside an "
    "analytic_geometry ... end block. The 'brick' keyword is used to begin a "
    "brick ... end block. At this time a brick can only be described by the "
    "coordinates of its oppposite points by their x,y,z coordinates. Example:\n"
    "\tanalytical_geometry\n"
    "\t  box 0.0 0.0 0.0  1.1 1.2 1.3 end # opposite corners: x1 y1 z1 x2 y2 "
    "z2\n"
    "\tend";
  }
};
using brick = keyword<brick_info, b,r,i,c,k>;

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
  static const char* name() { return
    "Interaction by exchange with the mean";
  }
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
  static const char* name() { return
    "Dirichlet";
  }
  static const char* help() { return
    "Material mix model, 'mix_dir', is short for Dirichlet. It is a material "
    "mix model that explicitly satisfies the unit-sum requirement for all "
    "statistical samples.";
  }
};
using mix_dir = keyword<mix_dir_info,  m,i,x,'_',d,i,r >;

// Keyword 'mix_gendir'
struct mix_gendir_info {
  static const char* name() { return
    "Generalized Dirichlet";
  }
  static const char* help() { return
    "Material mix model, 'mix_gendir', is short for generalized Dirichlet. It "
    "is a material mix model that explicitly satisfies the unit-sum "
    "requirement for all statistical samples.";
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

// Dirichlet and generalized Dirichlet parameters
using dir_B = keyword<undefined_info,  b >;
using dir_S = keyword<undefined_info,  S >;
using dir_kappa = keyword<undefined_info,  k,a,p,p,a >;
using gendir_C = keyword<undefined_info,  C >;

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

// Statistics output interval
using stai = keyword<undefined_info,  s,t,a,i >;

// PDF output interval
using pdfi = keyword<undefined_info,  p,d,f,i >;

// Glob output interval
using glbi = keyword<undefined_info,  g,l,b,i >;

// Random number generator
using rng = keyword<undefined_info, r,n,g >;

// Random number generator seed
using seed = keyword<undefined_info, s,e,e,d >;

// Statistics
using statistics = keyword<undefined_info,  s,t,a,t,i,s,t,i,c,s >;

} // kw::
} // quinoa::

#undef Keywords

#endif // QuinoaInputDeckKeywords_h
