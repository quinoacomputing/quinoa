//******************************************************************************
/*!
  \file      src/Control/RNGTestKeywords.h
  \author    J. Bakosi
  \date      Thu Sep 19 09:34:18 2013
  \copyright Copyright 2005-2012, Jozsef Bakosi, All rights reserved.
  \brief     Random number generator test suite keywords
  \details   Random number generator test suite keywords
*/
//******************************************************************************
#ifndef RNGTestKeywords_h
#define RNGTestKeywords_h

//! Signal to compiler that we are building a list of keywords. This is used by
//! the inline includes below to make sure they get included in the correct
//! namespace and not polluting the global one.
#define Keywords

namespace rngtest {
namespace grm {
namespace kw {

  using namespace pegtl::ascii;

  // Include base keywords recognized by all parsers
  #include <BaseKeywords.h>

  // Random number generator (RNG) test suite start block keyword
  using rngtest = pegtl::string< r,n,g,t,e,s,t >;

  // RNG test suite selector keyword
  using suite = pegtl::string< s,u,i,t,e >;

  // Various RNG test suite keywords
  using smallcrush = pegtl::string< s,m,a,l,l,c,r,u,s,h >;
  using crush = pegtl::string< c,r,u,s,h >;
  using bigcrush = pegtl::string< b,i,g,c,r,u,s,h >;

  // List of RNGs start block keyword
  using rngs = pegtl::string< r,n,g,s >;

  // MKL RNGs
  using mkl_mcg31 = pegtl::string< m,k,l,'_',m,c,g,'3','1' >;
  using mkl_r250 = pegtl::string< m,k,l,'_',r,'2','5','0' >;
  using mkl_mrg32k3a = pegtl::string< m,k,l,'_',m,r,g,'3','2',k,'3',a >;
  using mkl_mcg59 = pegtl::string< m,k,l,'_',m,c,g,'5','9' >;
  using mkl_wh = pegtl::string< m,k,l,'_',w,h >;
  using mkl_mt19937 = pegtl::string< m,k,l,'_',m,t,'1','9','9','3','7' >;
  using mkl_mt2203 = pegtl::string< m,k,l,'_',m,t,'2','2','0','3' >;
  using mkl_sfmt19937 = pegtl::string< m,k,l,'_',s,f,m,t,'1','9','9','3','7' >;
  using mkl_sobol = pegtl::string< m,k,l,'_',s,o,b,o,l >;
  using mkl_niederr = pegtl::string< m,k,l,'_',n,i,e,d,e,r,r >;
  using mkl_iabstract = pegtl::string< m,k,l,'_',i,a,b,s,t,r,a,c,t >;
  using mkl_dabstract = pegtl::string< m,k,l,'_',d,a,b,s,t,r,a,c,t >;
  using mkl_sabstract = pegtl::string< m,k,l,'_',s,a,b,s,t,r,a,c,t >;
  using mkl_nondeterm = pegtl::string< m,k,l,'_',n,o,n,d,e,t,e,r,m >;

} // kw::
} // grm::
} // rngtest::

#undef Keywords

#endif // RNGTestKeywords_h
