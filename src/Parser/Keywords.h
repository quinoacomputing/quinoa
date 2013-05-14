//******************************************************************************
/*!
  \file      src/Parser/Keywords.h
  \author    J. Bakosi
  \date      Mon 13 May 2013 08:58:48 PM MDT
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

  //****************************************************************************
  // Editing anything in this section should be accompanied by the corresponding
  // changes in FwdAssociate.h and BackAssociate.h in src/Control.

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
  using invpos = pegtl::string< i,n,v,p,o,s >;
  //   * Viscous model
  using vispos = pegtl::string< v,i,s,p,o,s >;

  // Select mass model:
  //   * Beta model
  using beta = pegtl::string< b,e,t,a >;

  // Select hydrodynamics model:
  //   * Simplified Langevin model
  using slm = pegtl::string< s,l,m >;
  //   * Generalized Langevin model
  using glm = pegtl::string< g,l,m >;

  // Select material mix model:
  //   * Interaction by exchange with the mean
  using iem = pegtl::string< i,e,m >;
  //   * Interaction by exchange with the conditional mean
  using iecm = pegtl::string< i,e,c,m >;
  //   * Dirichlet
  using dir = pegtl::string< d,i,r >;
  //   * generalized Dirichlet
  using gendir = pegtl::string< g,e,n,d,i,r >;
  //****************************************************************************

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

  // Dirichlet and generalized Dirichlet parameters
  using dir_B = pegtl::string< b >;
  using dir_S = pegtl::string< S >;
  using dir_kappa = pegtl::string< k,a,p,p,a >;
  using gendir_C = pegtl::string< C >;

  // Langevin model parameters
  using SLM_C0 = pegtl::string< C,'0' >;

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
