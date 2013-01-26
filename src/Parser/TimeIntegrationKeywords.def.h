//******************************************************************************
/*!
  \file      src/Parser/TimeIntegrationKeywords.def.h
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 07:41:38 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Time integration keywords
  \details   Keywords recognized in the time_integration - end block
*/
//******************************************************************************

namespace keyword {

  using nsteps = pegtl::string< n,s,t,e,p,s >;
  using term   = pegtl::string< t, e, r, m >;
  using dt     = pegtl::string< d,t >;

  struct time_integration :
         sor< nsteps, term, dt > {};

} // namespace keyword
