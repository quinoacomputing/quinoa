//******************************************************************************
/*!
  \file      src/Control/MKLRNGKeywords.h
  \author    J. Bakosi
  \date      Fri 22 Nov 2013 05:20:02 PM MST
  \copyright 2005-2014, Jozsef Bakosi.
  \brief     RNG keywords for Intel's MKL
  \details   Random number generator selector keywords for those generators
  available in Intel's Math Kernel Library's (MKL) Vector Statistical Library
  (VSL). The keywords are defined by specializing struct 'keyword', defined in
  Control/Keyword.h. Introducing a new keyword requires a more human readable
  (bust still short) name as well as a short, few-line, help-like description.
*/
//******************************************************************************
#ifndef Keywords
#error "MKLRNGKeywords.h should only be included within a *Keywords.h"
#endif

namespace tk {
namespace kw {

using namespace pegtl::ascii;

// Keyword 'mkl_mcg31'
struct mkl_mcg31_info {
  static const char* name() { return "Intel MKL MCG31"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_MCG31', a 31-bit multiplicative "
    "congruential random number generator, available in Intel's Math Kernel "
    "Library (MKL).";
  }
};
using mkl_mcg31 =
  keyword< mkl_mcg31_info, m,k,l,'_',m,c,g,'3','1' >;

// Keyword 'mkl_r250'
struct mkl_r250_info {
  static const char* name() { return "Intel MKL R250"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_R250', a generalized feedback "
    "shift register random number generator, available in Intel's Math Kernel "
    "Library (MKL).";
  }
};
using mkl_r250 =
  keyword< mkl_r250_info, m,k,l,'_',r,'2','5','0' >;

// Keyword 'mkl_mrg32k3a'
struct mkl_mrg32k3a_info {
  static const char* name() { return "Intel MKL MRG32K3A"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_MRG32K3A', a combined multiple "
    "recursive random number generator with two components of order 3, "
    "available in Intel's Math Kernel Library (MKL).";
  }
};
using mkl_mrg32k3a =
  keyword< mkl_mrg32k3a_info, m,k,l,'_',m,r,g,'3','2',k,'3',a >;

// Keyword 'mkl_mcg59'
struct mkl_mcg59_info {
  static const char* name() { return "Intel MKL MCG59"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_MCG59', a 59-bit multiplicative "
    "congruential random number generator, available in Intel's Math Kernel "
    "Library (MKL).";
  }
};
using mkl_mcg59 = keyword< mkl_mcg59_info, m,k,l,'_',m,c,g,'5','9' >;

// Keyword 'mkl_wh'
struct mkl_wh_info {
  static const char* name() { return "Intel MKL WH"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_WH', a set of 273 Wichmann-Hill "
    "combined multiplicative congruential random number generators, available "
    "in Intel's Math Kernel Library (MKL).";
  }
};
using mkl_wh = keyword< mkl_wh_info, m,k,l,'_',w,h >;

// Keyword 'mkl_mt19937'
struct mkl_mt19937_info {
  static const char* name() { return "Intel MKL MT19937"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_MT19937', a Mersenne Twister "
    "pseudorandom number generator, available in Intel's Math Kernel Library "
    "(MKL).";
  }
};
using mkl_mt19937 =
  keyword< mkl_mt19937_info, m,k,l,'_',m,t,'1','9','9','3','7' >;

// Keyword 'mkl_mt2203'
struct mkl_mt2203_info {
  static const char* name() { return "Intel MKL MT2203"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_MT2203', a set of 6024 Mersenne "
    "Twister pseudorandom number generators, available in Intel's Math Kernel "
    "Library (MKL).";
  }
};
using mkl_mt2203 = keyword< mkl_mt2203_info, m,k,l,'_',m,t,'2','2','0','3' >;

// Keyword 'mkl_sfmt19937'
struct mkl_sfmt19937_info {
  static const char* name() { return "Intel MKL SFMT19937"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_SFMT19937', a SIMD-oriented Fast "
    "Mersenne Twister pseudorandom number generator, available in Intel's Math "
    "Kernel Library (MKL).";
  }
};
using mkl_sfmt19937 =
  keyword< mkl_sfmt19937_info, m,k,l,'_',s,f,m,t,'1','9','9','3','7' >;

// Keyword 'mkl_sobol'
struct mkl_sobol_info {
  static const char* name() { return "Intel MKL SOBOL"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_SOBOL', a 32-bit Gray code-based "
    "random number generator, producing low-discrepancy sequences for "
    "dimensions 1 ≤ s ≤ 40 with available user-defined dimensions, available "
    "in Intel's Math Kernel Library (MKL).";
  }
};
using mkl_sobol = keyword< mkl_sobol_info, m,k,l,'_',s,o,b,o,l >;

// Keyword 'mkl_niederr'
struct mkl_niederr_info {
  static const char* name() { return "Intel MKL NIEDERR"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_NIEDERR', a 32-bit Gray "
    "code-based random number generator, producing low-discrepancy sequences "
    "for dimensions 1 ≤ s ≤ 318 with available user-defined dimensions, "
    "available in Intel's Math Kernel Library (MKL).";
  }
};
using mkl_niederr = keyword< mkl_niederr_info, m,k,l,'_',n,i,e,d,e,r,r >;

// Keyword 'mkl_iabstract'
struct mkl_iabstract_info {
  static const char* name() { return "Intel MKL IABSTRACT"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_IABSTRACT', an abstract random "
    "number generator for integer arrays, available in Intel's Math Kernel "
    "Library (MKL).";
  }
};
using mkl_iabstract =
  keyword< mkl_iabstract_info, m,k,l,'_',i,a,b,s,t,r,a,c,t >;

// Keyword 'mkl_dabstract'
struct mkl_dabstract_info {
  static const char* name() { return "Intel MKL DABSTRACT"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_DABSTRACT', an abstract random "
    "number generator for double-precision floating-point arrays, available in "
    "Intel's Math Kernel Library (MKL).";
  }
};
using mkl_dabstract =
  keyword< mkl_dabstract_info, m,k,l,'_',d,a,b,s,t,r,a,c,t >;

// Keyword 'mkl_sabstract'
struct mkl_sabstract_info {
  static const char* name() { return "Intel MKL SABSTRACT"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_SABSTRACT', an abstract random "
    "number generator for single-precision floating-point arrays, available in "
    "Intel's Math Kernel Library (MKL).";
  }
};
using mkl_sabstract =
  keyword< mkl_sabstract_info, m,k,l,'_',s,a,b,s,t,r,a,c,t >;

// Keyword 'mkl_nondeterm'
struct mkl_nondeterm_info {
  static const char* name() { return "Intel MKL NONDETERM"; }
  static const char* help() { return
    "This keyword is used to select 'VSL_BRNG_NONDETERM', a non-deterministic "
    "random number generator, available in Intel's Math Kernel Library (MKL).";
  }
};
using mkl_nondeterm =
  keyword< mkl_nondeterm_info, m,k,l,'_',n,o,n,d,e,t,e,r,m >;

// Keyword 'uniform_method'
struct uniform_method_info {
  static const char* name() { return "Intel MKL uniform RNG method"; }
  static const char* help() { return
    "This keyword is used to specify the method used to generate uniform "
    "random numbers using the Intel MKL library. Valid options are 'std' and "
    "'accurate'.";
  }
};
using uniform_method = keyword< uniform_method_info,
                                u,n,i,f,o,r,m,'_',m,e,t,h,o,d >;

// Keyword 'standard'
struct standard_info {
  static const char* name() { return "standard"; }
  static const char* help() { return
    "This keyword is used to specify the standard method used to generate "
    "uniform random numbers using the Intel MKL library.";
  }
};
using standard = keyword< standard_info, s,t,a,n,d,a,r,d >;

// Keyword 'accurate'
struct accurate_info {
  static const char* name() { return "accurate"; }
  static const char* help() { return
    "This keyword is used to specify the accurate method used to generate "
    "uniform random numbers using the Intel MKL library.";
  }
};
using accurate = keyword< accurate_info, a,c,c,u,r,a,t,e >;

// Keyword 'gaussian_method'
struct gaussian_method_info {
  static const char* name() { return "Intel MKL Gaussian RNG method"; }
  static const char* help() { return
    "This keyword is used to specify the method used to generate Gaussian "
    "random numbers using the Intel MKL library. Valid options are 'boxmuller' "
    ", 'boxmuller2', and 'icdf'.";
  }
};
using gaussian_method = keyword< gaussian_method_info,
                                 g,a,u,s,s,i,a,n,'_',m,e,t,h,o,d >;

// Keyword 'boxmuller'
struct boxmuller_info {
  static const char* name() { return "boxmuller"; }
  static const char* help() { return
    "This keyword is used to specify the Box-Muller method used to generate "
    "Gaussian random numbers using the Intel MKL library.";
  }
};
using boxmuller = keyword< boxmuller_info, b,o,x,m,u,l,l,e,r >;

// Keyword 'boxmuller2'
struct boxmuller2_info {
  static const char* name() { return "boxmuller2"; }
  static const char* help() { return
    "This keyword is used to specify the Box-Muller 2 method used to generate "
    "Gaussian random numbers using the Intel MKL library.";
  }
};
using boxmuller2 = keyword< boxmuller2_info, b,o,x,m,u,l,l,e,r,'2' >;

// Keyword 'icdf'
struct icdf_info {
  static const char* name() { return "icdf"; }
  static const char* help() { return
    "This keyword is used to specify the ICDF method used to generate "
    "Gaussian random numbers using the Intel MKL library.";
  }
};
using icdf = keyword< icdf_info, i,c,d,f >;

} // kw::
} // tk::
