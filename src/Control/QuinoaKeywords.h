//******************************************************************************
/*!
  \file      src/Control/QuinoaKeywords.h
  \author    J. Bakosi
  \date      Sat 21 Sep 2013 07:37:40 AM MDT
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Quinoa's keywords
  \details   All keywords recognized by Quinoa's parser
*/
//******************************************************************************
#ifndef QuinoaKeywords_h
#define QuinoaKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes below to make sure they get included in the correct
//! namespace and not polluting the global one.
#define Keywords

#include <pegtl.hh>

namespace quinoa {
namespace grm {

//! List of keywords the parser understands
namespace kw {

//! A keyword is struct that combines a type (pegtl::string) and a value
//! (std::string)
template< const char* Name, const char* Help, int Char, int... Chars >
struct keyword {

  //! Accessor as pegtl::tring
  using pegtl_string = pegtl::string<Char, Chars...>;

  //! Accessor as std::tring
  std::string string() const {
    return std::move(std::string( (sizeof...(Chars)) ?
                                  (pegtl::escaper<Char, Chars...>::result()) :
                                  (pegtl::escape(Char)) ));
  }

  //! Accessor to name
  const char* name() const { return Name; }
  //! Accessor to help
  const char* help() const { return Help; }
};

using namespace pegtl::ascii;

// Include base keywords recognized by all parsers
#include <BaseKeywords.h>

// Keyword 'analytic_geometry'
const char ag_name[] = "Analytic geometry defintion";
const char ag_help[] =
  "This option is used to define the beginning of an analytical geometry "
  "definition block. Analytical geometry definitions describe a 3D geometry "
  "using available primitives. Example:\n"
  "\tanalytical_geometry\n"
  "\t  box 0.0 0.0 0.0  1.1 1.2 1.3 end # opposite corners: x1 y1 z1 x2 y2 z2\n"
  "\tend";
using analytic_geometry =
  keyword<ag_name, ag_help, a,n,a,l,y,t,i,c,'_',g,e,o,m,e,t,r,y>;

// Keyword 'discrete_geometry'
const char dg_name[] = "Discrete geometry defintion";
const char dg_help[] =
  "This option is used to define the beginning of a discrete geometry "
  "definition block. Discrete geometry definitions are used to specify various "
  "parameters that influence the generation of point clouds within a closed "
  "geometry specified by discrete surfaces. Example:\n"
  "\tdiscrete_geometry\n"
  "\t  ...\n"
  "\tend";
using discrete_geometry =
  keyword<ag_name, ag_help, d,i,s,c,r,e,t,e,'_',g,e,o,m,e,t,r,y>;

// Geometry primitives for analytic geometry definition
//   * Box
using box = keyword<ag_name, ag_help, b,o,x>;

// Input filename
using input = keyword<ag_name, ag_help, i,n,p,u,t>;

// Output filename
using output = keyword<ag_name, ag_help, o,u,t,p,u,t>;

// PDF filename
using pdfname = keyword<ag_name, ag_help,  p,d,f,n,a,m,e >;

// Glob (i.e. domain-average statistics) filename
using globname = keyword<ag_name, ag_help,  g,l,o,b,n,a,m,e >;

// Statistics filename
using statname = keyword<ag_name, ag_help,  s,t,a,t,n,a,m,e >;

// Select physics:
//   * Homogeneous material mixing
using hommix = keyword<ag_name, ag_help, h,o,m,m,i,x>;
//   * Homogeneous hydrodinamics
using homhydro = keyword<ag_name, ag_help,  h,o,m,h,y,d,r,o >;
//   * Homogeneous Rayleigh-Taylor
using homrt = keyword<ag_name, ag_help,  h,o,m,r,t >;
//   * Standalone-particle incompressible Navier-Stokes flow
using spinsflow = keyword<ag_name, ag_help,  s,p,i,n,s,f,l,o,w >;

// Select position model:
//   * Insviscid model
using pos_inviscid = keyword<ag_name, ag_help,  p,o,s,'_',i,n,v,i,s,c,i,d >;
//   * Viscous model
using pos_viscous = keyword<ag_name, ag_help,  p,o,s,'_',v,i,s,c,o,u,s >;

// Select mass model:
//   * Beta model
using mass_beta = keyword<ag_name, ag_help,  m,a,s,s,'_',b,e,t,a >;

// Select hydrodynamics model:
//   * Simplified Langevin model
using hydro_slm = keyword<ag_name, ag_help,  h,y,d,r,o,'_',s,l,m >;
//   * Generalized Langevin model
using hydro_glm = keyword<ag_name, ag_help,  h,y,d,r,o,'_',g,l,m >;

// Select material mix model:
//   * Interaction by exchange with the mean
using mix_iem = keyword<ag_name, ag_help,  m,i,x,'_',i,e,m >;
//   * Interaction by exchange with the conditional mean
using mix_iecm = keyword<ag_name, ag_help,  m,i,x,'_',i,e,c,m >;
//   * Dirichlet
using mix_dir = keyword<ag_name, ag_help,  m,i,x,'_',d,i,r >;
//   * generalized Dirichlet
using mix_gendir = keyword<ag_name, ag_help,  m,i,x,'_',g,e,n,d,i,r >;

// Select material mix rate model:
//   * Gamma distribution model
using mixrate_gamma = keyword<ag_name, ag_help,  m,i,x,r,a,t,e,'_',g,a,m,m,a >;

// Select turbulence frequency model:
//   * Gamma distribution model
using freq_gamma = keyword<ag_name, ag_help,  f,r,e,q,'_',g,a,m,m,a >;

// Number of time steps to take
using nstep = keyword<ag_name, ag_help,  n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword<ag_name, ag_help,  t, e, r, m >;

// Size of time step
using dt = keyword<ag_name, ag_help,  d,t >;

// Start of position model specification block
using position = keyword<ag_name, ag_help,  p,o,s,i,t,i,o,n >;

// Start of hydrodynamics model specification block
using hydro = keyword<ag_name, ag_help,  h,y,d,r,o >;

// Start of material mix model specification block
using mix = keyword<ag_name, ag_help,  m,i,x >;

// Number of particle position components
using nposition = keyword<ag_name, ag_help,  n,p,o,s,i,t,i,o,n >;
// Number of particle density components
using ndensity = keyword<ag_name, ag_help,  n,d,e,n,s,i,t,y >;
// Number of particle velocity components
using nvelocity = keyword<ag_name, ag_help,  n,v,e,l,o,c,i,t,y >;
// Number of particle scalar components
using nscalar = keyword<ag_name, ag_help,  n,s,c,a,l,a,r >;
// Number of particle turbulence frequency components
using nfreq = keyword<ag_name, ag_help,  n,f,r,e,q >;

// Dirichlet and generalized Dirichlet parameters
using dir_B = keyword<ag_name, ag_help,  b >;
using dir_S = keyword<ag_name, ag_help,  S >;
using dir_kappa = keyword<ag_name, ag_help,  k,a,p,p,a >;
using gendir_C = keyword<ag_name, ag_help,  C >;

// Langevin model parameters
using SLM_C0 = keyword<ag_name, ag_help,  C,'0' >;

// Gamma frequency model parameters
using freq_gamma_C1 = keyword<ag_name, ag_help,  C,'1' >;
using freq_gamma_C2 = keyword<ag_name, ag_help,  C,'2' >;
using freq_gamma_C3 = keyword<ag_name, ag_help,  C,'3' >;
using freq_gamma_C4 = keyword<ag_name, ag_help,  C,'4' >;

// Beta model parameters
using Beta_At = keyword<ag_name, ag_help,  A,t >;

// Quantities
using transported_scalar = keyword<ag_name, ag_help,  Y >;
using transported_scalar_fluctuation = keyword<ag_name, ag_help,  y >;

using velocity_x = keyword<ag_name, ag_help,  U >;
using velocity_fluctuation_x = keyword<ag_name, ag_help,  u >;
using velocity_y = keyword<ag_name, ag_help,  V >;
using velocity_fluctuation_y = keyword<ag_name, ag_help,  v >;
using velocity_z = keyword<ag_name, ag_help,  W >;
using velocity_fluctuation_z = keyword<ag_name, ag_help,  w >;

using pressure = keyword<ag_name, ag_help,  P >;
using pressure_fluctuation = keyword<ag_name, ag_help,  p >;
  
using density = keyword<ag_name, ag_help,  R >;
using density_fluctuation = keyword<ag_name, ag_help,  r >;
  
// Total number of particles
using npar = keyword<ag_name, ag_help,  n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword<ag_name, ag_help,  t,t,y,i >;

// Dump (restart file) output interval
using dmpi = keyword<ag_name, ag_help,  d,m,p,i >;

// Statistics output interval
using stai = keyword<ag_name, ag_help,  s,t,a,i >;

// PDF output interval
using pdfi = keyword<ag_name, ag_help,  p,d,f,i >;

// Glob output interval
using glbi = keyword<ag_name, ag_help,  g,l,b,i >;

// Statistics
using statistics = keyword<ag_name, ag_help,  s,t,a,t,i,s,t,i,c,s >;

// Random number generator (RNG) test suite
using rngtest = keyword<ag_name, ag_help,  r,n,g,t,e,s,t >;

// RNG test suite
using suite = keyword<ag_name, ag_help,  s,u,i,t,e >;

// RNG test suites
using smallcrush = keyword<ag_name, ag_help,  s,m,a,l,l,c,r,u,s,h >;
using crush = keyword<ag_name, ag_help,  c,r,u,s,h >;
using bigcrush = keyword<ag_name, ag_help,  b,i,g,c,r,u,s,h >;

// RNGs
using rngs = keyword<ag_name, ag_help,  r,n,g,s >;

// MKL RNGs
using mkl_mcg31 = keyword<ag_name, ag_help,  m,k,l,'_',m,c,g,'3','1' >;
using mkl_r250 = keyword<ag_name, ag_help,  m,k,l,'_',r,'2','5','0' >;
using mkl_mrg32k3a = keyword<ag_name, ag_help,  m,k,l,'_',m,r,g,'3','2',k,'3',a >;
using mkl_mcg59 = keyword<ag_name, ag_help,  m,k,l,'_',m,c,g,'5','9' >;
using mkl_wh = keyword<ag_name, ag_help,  m,k,l,'_',w,h >;
using mkl_mt19937 = keyword<ag_name, ag_help,  m,k,l,'_',m,t,'1','9','9','3','7' >;
using mkl_mt2203 = keyword<ag_name, ag_help,  m,k,l,'_',m,t,'2','2','0','3' >;
using mkl_sfmt19937 = keyword<ag_name, ag_help,  m,k,l,'_',s,f,m,t,'1','9','9','3','7' >;
using mkl_sobol = keyword<ag_name, ag_help,  m,k,l,'_',s,o,b,o,l >;
using mkl_niederr = keyword<ag_name, ag_help,  m,k,l,'_',n,i,e,d,e,r,r >;
using mkl_iabstract = keyword<ag_name, ag_help,  m,k,l,'_',i,a,b,s,t,r,a,c,t >;
using mkl_dabstract = keyword<ag_name, ag_help,  m,k,l,'_',d,a,b,s,t,r,a,c,t >;
using mkl_sabstract = keyword<ag_name, ag_help,  m,k,l,'_',s,a,b,s,t,r,a,c,t >;
using mkl_nondeterm = keyword<ag_name, ag_help,  m,k,l,'_',n,o,n,d,e,t,e,r,m >;

} // kw::
} // grm::
} // quinoa::

#endif // QuinoaKeywords_h
