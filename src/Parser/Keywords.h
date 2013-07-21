//******************************************************************************
/*!
  \file      src/Parser/Keywords.h
  \author    J. Bakosi
  \date      Fri Jul 19 16:05:03 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Keywords
  \details   All keywords recognized by the parser
*/
//******************************************************************************
#ifndef Grammar_h
#error "Keywords.h should only be included within Grammar.h"
#endif

#ifndef Keywords_h
#define Keywords_h

// Keywords accepted by the parser
namespace keyword {

  // Problem title
  using title = pegtl::string<t,i,t,l,e>;

  // End of block
  using end = pegtl::string< e,n,d >;

  // Select geometry definition
  //   * Analytic
  using analytic_geometry = pegtl::string<a,n,a,l,y,t,i,c,'_',g,e,o,m,e,t,r,y>;
  //   * Discrete
  using discrete_geometry = pegtl::string<d,i,s,c,r,e,t,e,'_',g,e,o,m,e,t,r,y>;

  // Geometry primitives for analytic geometry definition
  //   * Box
  using box = pegtl::string<b,o,x>;

  // Geometry input filename
  using input = pegtl::string<i,n,p,u,t>;
  // Geometry output filename
  using output = pegtl::string<o,u,t,p,u,t>;

  // Select physics:
  //   * Homogeneous material mixing
  using hommix = pegtl::string< h,o,m,m,i,x >;
  //   * Homogeneous hydrodinamics
  using homhydro = pegtl::string< h,o,m,h,y,d,r,o >;
  //   * Homogeneous Rayleigh-Taylor
  using homrt = pegtl::string< h,o,m,r,t >;
  //   * Standalone-particle incompressible Navier-Stokes flow
  using spinsflow = pegtl::string< s,p,i,n,s,f,l,o,w >;

  // Select position model:
  //   * Insviscid model
  using pos_inviscid = pegtl::string< p,o,s,'_',i,n,v,i,s,c,i,d >;
  //   * Viscous model
  using pos_viscous = pegtl::string< p,o,s,'_',v,i,s,c,o,u,s >;

  // Select mass model:
  //   * Beta model
  using mass_beta = pegtl::string< m,a,s,s,'_',b,e,t,a >;

  // Select hydrodynamics model:
  //   * Simplified Langevin model
  using hydro_slm = pegtl::string< h,y,d,r,o,'_',s,l,m >;
  //   * Generalized Langevin model
  using hydro_glm = pegtl::string< h,y,d,r,o,'_',g,l,m >;

  // Select material mix model:
  //   * Interaction by exchange with the mean
  using mix_iem = pegtl::string< m,i,x,'_',i,e,m >;
  //   * Interaction by exchange with the conditional mean
  using mix_iecm = pegtl::string< m,i,x,'_',i,e,c,m >;
  //   * Dirichlet
  using mix_dir = pegtl::string< m,i,x,'_',d,i,r >;
  //   * generalized Dirichlet
  using mix_gendir = pegtl::string< m,i,x,'_',g,e,n,d,i,r >;

  // Select turbulence frequency model:
  //   * Gamma distribution model
  using freq_gamma = pegtl::string< f,r,e,q,'_',g,a,m,m,a >;

  // Number of time steps to take
  using nstep = pegtl::string< n,s,t,e,p >;

  // Terminate time stepping at this value
  using term = pegtl::string< t, e, r, m >;

  // Size of time step
  using dt = pegtl::string< d,t >;

  // Start of position model specification block
  using position = pegtl::string< p,o,s,i,t,i,o,n >;

  // Start of hydrodynamics model specification block
  using hydro = pegtl::string< h,y,d,r,o >;

  // Start of material mix model specification block
  using mix = pegtl::string< m,i,x >;

  // Number of particle position components
  using nposition = pegtl::string< n,p,o,s,i,t,i,o,n >;
  // Number of particle density components
  using ndensity = pegtl::string< n,d,e,n,s,i,t,y >;
  // Number of particle velocity components
  using nvelocity = pegtl::string< n,v,e,l,o,c,i,t,y >;
  // Number of particle scalar components
  using nscalar = pegtl::string< n,s,c,a,l,a,r >;
  // Number of particle turbulence frequency components
  using nfreq = pegtl::string< n,f,r,e,q >;

  // Dirichlet and generalized Dirichlet parameters
  using dir_B = pegtl::string< b >;
  using dir_S = pegtl::string< S >;
  using dir_kappa = pegtl::string< k,a,p,p,a >;
  using gendir_C = pegtl::string< C >;

  // Langevin model parameters
  using SLM_C0 = pegtl::string< C,'0' >;

  // Gamma frequency model parameters
  using freq_gamma_C1 = pegtl::string< C,'1' >;
  using freq_gamma_C2 = pegtl::string< C,'2' >;
  using freq_gamma_C3 = pegtl::string< C,'3' >;
  using freq_gamma_C4 = pegtl::string< C,'4' >;

  // Beta model parameters
  using Beta_At = pegtl::string< A,t >;

  // Quantities
  using transported_scalar = pegtl::string< Y >;
  using transported_scalar_fluctuation = pegtl::string< y >;

  using velocity_x = pegtl::string< U >;
  using velocity_fluctuation_x = pegtl::string< u >;
  using velocity_y = pegtl::string< V >;
  using velocity_fluctuation_y = pegtl::string< v >;
  using velocity_z = pegtl::string< W >;
  using velocity_fluctuation_z = pegtl::string< w >;

  using pressure = pegtl::string< P >;
  using pressure_fluctuation = pegtl::string< p >;
  
  using density = pegtl::string< R >;
  using density_fluctuation = pegtl::string< r >;
  
  // Total number of particles
  using npar = pegtl::string< n,p,a,r >;

  // TTY (screen) output interval
  using ttyi = pegtl::string< t,t,y,i >;

  // Dump (restart file) output interval
  using dump = pegtl::string< d,u,m,p >;

  // Plot output interval
  using plti = pegtl::string< p,l,t,i >;

  // PDF output interval
  using pdfi = pegtl::string< p,d,f,i >;

  // Glob output interval
  using glob = pegtl::string< g,l,o,b >;

  // Statistics
  using statistics = pegtl::string< s,t,a,t,i,s,t,i,c,s >;

  // PDF base filename
  using pdfname = pegtl::string< p,d,f,n,a,m,e >;

  // Glob (i.e. domain-average statistics) filename
  using globname = pegtl::string< g,l,o,b,n,a,m,e >;

  // Plot base filename
  using plotname = pegtl::string< p,l,o,t,n,a,m,e >;

} // namespace keyword

#endif // Keywords_h
