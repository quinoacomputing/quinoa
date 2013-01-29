//******************************************************************************
/*!
  \file      src/Parser/Keywords.def.h
  \author    J. Bakosi
  \date      Tue 29 Jan 2013 06:40:26 AM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Keywords
  \details   All keywords recognized by the parser
*/
//******************************************************************************

namespace keyword {

  // Main keywords
  using title = pegtl::string<t,i,t,l,e>;
  using end   = pegtl::string< e,n,d >;

  // Physics keywords
  using homdir    = pegtl::string< h,o,m,d,i,r >;
  using homgendir = pegtl::string< h,o,m,g,e,n,d,i,r >;
  using spinsflow = pegtl::string< s,p,i,n,s,f,l,o,w >;

  // Hydro model keywords
  using slm = pegtl::string< s,l,m >;
  using glm = pegtl::string< g,l,m >;

  // Mix model keywords
  using iem    = pegtl::string< i,e,m >;
  using iecm   = pegtl::string< i,e,c,m >;
  using dir    = pegtl::string< d,i,r >;
  using gendir = pegtl::string< g,e,n,d,i,r >;

  // Time integration keywords
  using nsteps = pegtl::string< n,s,t,e,p,s >;
  using term   = pegtl::string< t, e, r, m >;
  using dt     = pegtl::string< d,t >;

  // Physics::SPINSFlow keywords
  using hydro = pegtl::string< h,y,d,r,o >;
  using mix   = pegtl::string< m,i,x >;

  // Physics::HomogeneousDirichlet keywords
  using nscalar = pegtl::string< n,s,c,a,l,a,r >;
  using npar = pegtl::string< n,p,a,r >;
  using term = pegtl::string< t,e,r,m >;
  using nstep = pegtl::string< n,s,t,e,p >;
  using echo = pegtl::string< e,c,h,o >;

} // namespace keyword
