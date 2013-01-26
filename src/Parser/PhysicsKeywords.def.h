//******************************************************************************
/*!
  \file      src/Parser/PhysicsKeywords.def.h
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 07:39:31 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Physics keywords
  \details   Keywords selecting physics
*/
//******************************************************************************

namespace keyword {

  using HomogeneousDirichlet = pegtl::string< H,o,m,o,g,e,n,e,o,u,s,
                                              D,i,r,i,c,h,l,e,t >;
  using HomogeneousGeneralizedDirichlet = pegtl::string< H,o,m,o,g,e,n,e,o,u,s,
                                                         G,e,n,e,r,a,l,i,z,e,d,
                                                         D,i,r,i,c,h,l,e,t >;
  using SPINSFlow = pegtl::string<S,P,I,N,S,F,l,o,w>;

  struct physics :
         sor< HomogeneousDirichlet,
              HomogeneousGeneralizedDirichlet,
              SPINSFlow > {};

} // namespace keyword
