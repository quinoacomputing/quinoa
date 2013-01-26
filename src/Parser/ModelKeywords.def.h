//******************************************************************************
/*!
  \file      src/Parser/ModelKeywords.def.h
  \author    J. Bakosi
  \date      Fri 25 Jan 2013 07:41:28 PM MST
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Model keywords
  \details   Keywords selecting models
*/
//******************************************************************************

namespace keyword {

  using namespace pegtl::ascii;

  // Hydro model keywords

  using SimplifiedLangevin = pegtl::string< S,i,m,p,l,i,f,i,e,d,
                                            L,a,n,g,e,v,i,n >;
  using GeneralizedLangevin = pegtl::string< G,e,n,e,r,a,l,i,z,e,d,
                                             L,a,n,g,e,v,i,n >;

  struct hydro_model :
         sor< SimplifiedLangevin, GeneralizedLangevin > {};

  // Mix model keywords

  using Dirichlet = pegtl::string< D,i,r,i,c,h,l,e,t >;
  using GeneralizedDirichlet = pegtl::string< G,e,n,e,r,a,l,i,z,e,d,
                                              D,i,r,i,c,h,l,e,t >;

  struct mix_model :
         sor< Dirichlet, GeneralizedDirichlet > {};

  struct model :
         sor< hydro_model, mix_model > {};

} // namespace keyword
