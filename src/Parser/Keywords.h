//******************************************************************************
/*!
  \file      src/Parser/Keywords.h
  \author    J. Bakosi
  \date      Tue 29 Jan 2013 08:14:53 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Keywords
  \details   All keywords recognized by the parser
*/
//******************************************************************************

namespace keyword {

  // Problem title
  using title = pegtl::string<t,i,t,l,e>;

  // End of block
  using end = pegtl::string< e,n,d >;

  // Select physics:
  //   * Homogeneous Dirichlet
  using homdir = pegtl::string< h,o,m,d,i,r >;
  //   * Homogeneous generalized Dirichlet
  using homgendir = pegtl::string< h,o,m,g,e,n,d,i,r >;
  //   * Standalone-particle incompressible Navier-Stokes flow
  using spinsflow = pegtl::string< s,p,i,n,s,f,l,o,w >;

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

  // Number of time steps to take
  using nstep = pegtl::string< n,s,t,e,p >;

  // Terminate time stepping at this value
  using term = pegtl::string< t, e, r, m >;

  // Size of time step
  using dt = pegtl::string< d,t >;

  // Start of hydrodynamics model specification block
  using hydro = pegtl::string< h,y,d,r,o >;

  // Start of material mix model specification block
  using mix = pegtl::string< m,i,x >;

  // Number of mixing scalars
  using nscalar = pegtl::string< n,s,c,a,l,a,r >;

  // Total number of particles
  using npar = pegtl::string< n,p,a,r >;

  // One-line screen output every few time step
  using echo = pegtl::string< e,c,h,o >;

} // namespace keyword
