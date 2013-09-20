//******************************************************************************
/*!
  \file      src/Control/QuinoaKeywords.h
  \author    J. Bakosi
  \date      Fri Sep 20 13:40:14 2013
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
template< int Char, int... Chars >
struct keyword {
  using pegtl_string = pegtl::string<Char, Chars...>;
  std::string string() const {
    return std::move(std::string( (sizeof...(Chars)) ?
                        (pegtl::escaper<Char, Chars...>::result()) :
                        (pegtl::escape(Char)) ));
  }
};

using namespace pegtl::ascii;

// Include base keywords recognized by all parsers
#include <BaseKeywords.h>

// Select geometry definition
//   * Analytic
using analytic_geometry = keyword<a,n,a,l,y,t,i,c,'_',g,e,o,m,e,t,r,y>;
//   * Discrete
using discrete_geometry = keyword<d,i,s,c,r,e,t,e,'_',g,e,o,m,e,t,r,y>;

// Geometry primitives for analytic geometry definition
//   * Box
using box = keyword<b,o,x>;

// Input filename
using input = keyword<i,n,p,u,t>;

// Output filename
using output = keyword<o,u,t,p,u,t>;

// PDF filename
using pdfname = keyword< p,d,f,n,a,m,e >;

// Glob (i.e. domain-average statistics) filename
using globname = keyword< g,l,o,b,n,a,m,e >;

// Statistics filename
using statname = keyword< s,t,a,t,n,a,m,e >;

// Select physics:
//   * Homogeneous material mixing
using hommix = keyword<h,o,m,m,i,x>;
//   * Homogeneous hydrodinamics
using homhydro = keyword< h,o,m,h,y,d,r,o >;
//   * Homogeneous Rayleigh-Taylor
using homrt = keyword< h,o,m,r,t >;
//   * Standalone-particle incompressible Navier-Stokes flow
using spinsflow = keyword< s,p,i,n,s,f,l,o,w >;

// Select position model:
//   * Insviscid model
using pos_inviscid = keyword< p,o,s,'_',i,n,v,i,s,c,i,d >;
//   * Viscous model
using pos_viscous = keyword< p,o,s,'_',v,i,s,c,o,u,s >;

// Select mass model:
//   * Beta model
using mass_beta = keyword< m,a,s,s,'_',b,e,t,a >;

// Select hydrodynamics model:
//   * Simplified Langevin model
using hydro_slm = keyword< h,y,d,r,o,'_',s,l,m >;
//   * Generalized Langevin model
using hydro_glm = keyword< h,y,d,r,o,'_',g,l,m >;

// Select material mix model:
//   * Interaction by exchange with the mean
using mix_iem = keyword< m,i,x,'_',i,e,m >;
//   * Interaction by exchange with the conditional mean
using mix_iecm = keyword< m,i,x,'_',i,e,c,m >;
//   * Dirichlet
using mix_dir = keyword< m,i,x,'_',d,i,r >;
//   * generalized Dirichlet
using mix_gendir = keyword< m,i,x,'_',g,e,n,d,i,r >;

// Select material mix rate model:
//   * Gamma distribution model
using mixrate_gamma = keyword< m,i,x,r,a,t,e,'_',g,a,m,m,a >;

// Select turbulence frequency model:
//   * Gamma distribution model
using freq_gamma = keyword< f,r,e,q,'_',g,a,m,m,a >;

// Number of time steps to take
using nstep = keyword< n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword< t, e, r, m >;

// Size of time step
using dt = keyword< d,t >;

// Start of position model specification block
using position = keyword< p,o,s,i,t,i,o,n >;

// Start of hydrodynamics model specification block
using hydro = keyword< h,y,d,r,o >;

// Start of material mix model specification block
using mix = keyword< m,i,x >;

// Number of particle position components
using nposition = keyword< n,p,o,s,i,t,i,o,n >;
// Number of particle density components
using ndensity = keyword< n,d,e,n,s,i,t,y >;
// Number of particle velocity components
using nvelocity = keyword< n,v,e,l,o,c,i,t,y >;
// Number of particle scalar components
using nscalar = keyword< n,s,c,a,l,a,r >;
// Number of particle turbulence frequency components
using nfreq = keyword< n,f,r,e,q >;

// Dirichlet and generalized Dirichlet parameters
using dir_B = keyword< b >;
using dir_S = keyword< S >;
using dir_kappa = keyword< k,a,p,p,a >;
using gendir_C = keyword< C >;

// Langevin model parameters
using SLM_C0 = keyword< C,'0' >;

// Gamma frequency model parameters
using freq_gamma_C1 = keyword< C,'1' >;
using freq_gamma_C2 = keyword< C,'2' >;
using freq_gamma_C3 = keyword< C,'3' >;
using freq_gamma_C4 = keyword< C,'4' >;

// Beta model parameters
using Beta_At = keyword< A,t >;

// Quantities
using transported_scalar = keyword< Y >;
using transported_scalar_fluctuation = keyword< y >;

using velocity_x = keyword< U >;
using velocity_fluctuation_x = keyword< u >;
using velocity_y = keyword< V >;
using velocity_fluctuation_y = keyword< v >;
using velocity_z = keyword< W >;
using velocity_fluctuation_z = keyword< w >;

using pressure = keyword< P >;
using pressure_fluctuation = keyword< p >;
  
using density = keyword< R >;
using density_fluctuation = keyword< r >;
  
// Total number of particles
using npar = keyword< n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword< t,t,y,i >;

// Dump (restart file) output interval
using dmpi = keyword< d,m,p,i >;

// Statistics output interval
using stai = keyword< s,t,a,i >;

// PDF output interval
using pdfi = keyword< p,d,f,i >;

// Glob output interval
using glbi = keyword< g,l,b,i >;

// Statistics
using statistics = keyword< s,t,a,t,i,s,t,i,c,s >;

// Random number generator (RNG) test suite
using rngtest = keyword< r,n,g,t,e,s,t >;

// RNG test suite
using suite = keyword< s,u,i,t,e >;

// RNG test suites
using smallcrush = keyword< s,m,a,l,l,c,r,u,s,h >;
using crush = keyword< c,r,u,s,h >;
using bigcrush = keyword< b,i,g,c,r,u,s,h >;

// RNGs
using rngs = keyword< r,n,g,s >;

// MKL RNGs
using mkl_mcg31 = keyword< m,k,l,'_',m,c,g,'3','1' >;
using mkl_r250 = keyword< m,k,l,'_',r,'2','5','0' >;
using mkl_mrg32k3a = keyword< m,k,l,'_',m,r,g,'3','2',k,'3',a >;
using mkl_mcg59 = keyword< m,k,l,'_',m,c,g,'5','9' >;
using mkl_wh = keyword< m,k,l,'_',w,h >;
using mkl_mt19937 = keyword< m,k,l,'_',m,t,'1','9','9','3','7' >;
using mkl_mt2203 = keyword< m,k,l,'_',m,t,'2','2','0','3' >;
using mkl_sfmt19937 = keyword< m,k,l,'_',s,f,m,t,'1','9','9','3','7' >;
using mkl_sobol = keyword< m,k,l,'_',s,o,b,o,l >;
using mkl_niederr = keyword< m,k,l,'_',n,i,e,d,e,r,r >;
using mkl_iabstract = keyword< m,k,l,'_',i,a,b,s,t,r,a,c,t >;
using mkl_dabstract = keyword< m,k,l,'_',d,a,b,s,t,r,a,c,t >;
using mkl_sabstract = keyword< m,k,l,'_',s,a,b,s,t,r,a,c,t >;
using mkl_nondeterm = keyword< m,k,l,'_',n,o,n,d,e,t,e,r,m >;

} // kw::
} // grm::
} // quinoa::

#endif // QuinoaKeywords_h
