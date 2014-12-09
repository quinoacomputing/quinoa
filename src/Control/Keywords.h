//******************************************************************************
/*!
  \file      src/Control/Keywords.h
  \author    J. Bakosi
  \date      Tue 09 Dec 2014 08:58:41 AM MST
  \copyright 2012-2014, Jozsef Bakosi.
  \brief     All keywords shared among different executables
  \details   All keywords shared among different executables
*/
//******************************************************************************
#ifndef Keywords_h
#define Keywords_h

#include <Keyword.h>

namespace kw {

using namespace pegtl::ascii;

// Keyword 'title'
struct title_info {
  static const char* name() { return "Analysis title"; }
  static const char* help() { return
    "The analysis title may be specified in the input file using the 'title' "
    "keyword. The 'title' keyword must be followed by a doubly-quoted string "
    "specifying the analysis title. Example: title \"Example problem\".";
  }
};
using title = keyword< title_info, t,i,t,l,e >;

// Keyword 'end'
struct end_info {
  static const char* name() { return "End of block"; }
  static const char* help() { return
    "The end of a block is given by the 'end' keyword.";
  }
};
using end = keyword< end_info, e,n,d >;

// Keyword 'seed'
struct seed_info {
  static const char* name() { return "seed"; }
  static const char* help() { return
    "This keyword is used to denote the random number generator seed, if "
    "applicable. Example:\n"
    "\trngsse_gm55\n"
    "\t  seed 1234\n"
    "\tend";
  }
};
using seed = keyword< seed_info, s,e,e,d >;

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


// Keyword 'rngsse_gm19'
struct rngsse_gm19_info {
  static const char* name() { return "RNGSSE GM19"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM19 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm19 =
  keyword< rngsse_gm19_info, r,n,g,s,s,e,'_',g,m,'1','9' >;

// Keyword 'rngsse_gm29'
struct rngsse_gm29_info {
  static const char* name() { return "RNGSSE GM29"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM29 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm29 =
  keyword< rngsse_gm29_info, r,n,g,s,s,e,'_',g,m,'2','9' >;

// Keyword 'rngsse_gm31'
struct rngsse_gm31_info {
  static const char* name() { return "RNGSSE GM31"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM31 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm31 =
  keyword< rngsse_gm31_info, r,n,g,s,s,e,'_',g,m,'3','1' >;

// Keyword 'rngsse_gm55'
struct rngsse_gm55_info {
  static const char* name() { return "RNGSSE GM55"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM55 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm55 =
  keyword< rngsse_gm55_info, r,n,g,s,s,e,'_',g,m,'5','5' >;

// Keyword 'rngsse_gm61'
struct rngsse_gm61_info {
  static const char* name() { return "RNGSSE GM61"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GM61 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gm61 =
  keyword< rngsse_gm61_info, r,n,g,s,s,e,'_',g,m,'6','1' >;

// Keyword 'rngsse_gq58.1'
struct rngsse_gq581_info {
  static const char* name() { return "RNGSSE GQ58.1"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.1 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq581 =
  keyword< rngsse_gq581_info, r,n,g,s,s,e,'_',g,q,'5','8','.','1' >;

// Keyword 'rngsse_gq58.3'
struct rngsse_gq583_info {
  static const char* name() { return "RNGSSE GQ58.3"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.3 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq583 =
  keyword< rngsse_gq583_info, r,n,g,s,s,e,'_',g,q,'5','8','.','3' >;

// Keyword 'rngsse_gq58.4'
struct rngsse_gq584_info {
  static const char* name() { return "RNGSSE GQ58.4"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's GQ58.4 generator, using a "
    "method based on parallel evolution of an ensemble of transformations of "
    "a two-dimensional torus.";
  }
};

using rngsse_gq584 =
  keyword< rngsse_gq584_info, r,n,g,s,s,e,'_',g,q,'5','8','.','4' >;

// Keyword 'rngsse_mt19937'
struct rngsse_mt19937_info {
  static const char* name() { return "RNGSSE MT19937"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's MT19937 generator, a "
    "Mersenne Twister generator.";
  }
};

using rngsse_mt19937 =
  keyword< rngsse_mt19937_info, r,n,g,s,s,e,'_',m,t,'1','9','9','3','7' >;

// Keyword 'rngsse_lfsr113'
struct rngsse_lfsr113_info {
  static const char* name() { return "RNGSSE LFSR113"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's LFSR113 generator.";
  }
};

using rngsse_lfsr113 =
  keyword< rngsse_lfsr113_info, r,n,g,s,s,e,'_',l,f,s,r,'1','1','3' >;

// Keyword 'rngsse_mrg32k3a'
struct rngsse_mrg32k3a_info {
  static const char* name() { return "RNGSSE MRG32K3A"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's MRG32K3A generator, a "
    "combined multiple recursive random number generator with two components "
    "of order 3.";
  }
};

using rngsse_mrg32k3a =
  keyword< rngsse_mrg32k3a_info, r,n,g,s,s,e,'_',m,r,g,'3','2',k,'3',a >;

// Keyword 'seqlen'
struct seqlen_info {
  static const char* name() { return "RNGSSE sequence length"; }
  static const char* help() { return
    "This keyword is used to select RNGSSE library's generator sequence "
    "length. Valid options are 'short', 'medium', and 'long'.";
  }
};

using seqlen =
  keyword< seqlen_info, s,e,q,l,e,n >;

// Keyword 'short'
struct seq_short_info {
  static const char* name() { return "short"; }
  static const char* help() { return
    "This keyword is used to select the short sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_short =
  keyword< seq_short_info, s,h,o,r,t >;

// Keyword 'medium'
struct seq_med_info {
  static const char* name() { return "medium"; }
  static const char* help() { return
    "This keyword is used to select the medium sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_med =
  keyword< seq_med_info, m,e,d,i,u,m >;

// Keyword 'long'
struct seq_long_info {
  static const char* name() { return "long"; }
  static const char* help() { return
    "This keyword is used to select the long sequence length for RNGSSE "
    "library's generator sequence.";
  }
};

using seq_long =
  keyword< seq_long_info, l,o,n,g >;


// pdfs - end block keyword
using pdfs = keyword<undefined_info, p,d,f,s >;

// Keyword 'filetype': selected PDF output file type
struct filetype_info {
  static const char* name() { return "PDF output file type"; }
  static const char* help() { return "PDF output file type."; }
};
using pdf_filetype = keyword< filetype_info, f,i,l,e,t,y,p,e >;

// Keyword 'txt': PDF output file type: txt
struct txt_info {
  static const char* name() { return "text"; }
  static const char* help() { return "PDF text output file."; }
};
using txt = keyword< txt_info, t,x,t >;

// Keyword 'gmshtxt': PDF output file type: text gmsh
struct gmshtxt_info {
  static const char* name() { return "gmshtxt"; }
  static const char* help() { return "PDF GMSH txt file output."; }
};
using gmshtxt = keyword< gmshtxt_info, g,m,s,h,t,x,t >;

// Keyword 'gmshbin': PDF output file type: binary gmsh
struct gmshbin_info {
  static const char* name() { return "gmshbin"; }
  static const char* help() { return "PDF GMSH binary file output."; }
};
using gmshbin = keyword< gmshbin_info, g,m,s,h,b,i,n >;

// Keyword 'exodusii': PDF output file type: Exodus II
struct exodusii_info {
  static const char* name() { return "exodusii"; }
  static const char* help() { return "PDF Exodus II binary file output."; }
};
using exodusii = keyword< exodusii_info, e,x,o,d,u,s,i,i >;

// Keyword 'policy': selected PDF output file policy
struct policy_info {
  static const char* name() { return "PDF output file policy"; }
  static const char* help() { return "PDF output file policy."; }
};
using pdf_policy = keyword< policy_info, p,o,l,i,c,y >;

// Keyword 'overwrite': PDF output file policy: overwrite
struct overwrite_info {
  static const char* name() { return "overwrite"; }
  static const char* help() { return "PDF output overwrites the same file "
  "containing a single time step."; }
};
using overwrite = keyword< overwrite_info, o,v,e,r,w,r,i,t,e >;

// Keyword 'multiple': PDF output file policy: multiple
struct multiple_info {
  static const char* name() { return "multiple"; }
  static const char* help() { return "PDF output creates new file for each "
  "time step."; }
};
using multiple = keyword< multiple_info, m,u,l,t,i,p,l,e >;

// Keyword 'evolution': PDF output file policy: evolution
struct evolution_info {
  static const char* name() { return "evolution"; }
  static const char* help() { return "PDF output appends new time step to the "
  "same file for each yielding a time evolution of data."; }
};
using evolution = keyword< evolution_info, e,v,o,l,u,t,i,o,n >;

// Keyword 'txt_float_format': PDF txt output floating-point format
struct txt_float_format_info {
  static const char* name() { return "format"; }
  static const char* help() { return "Formatting for floating-point output."; }
};
using txt_float_format = keyword< txt_float_format_info, f,o,r,m,a,t >;

// Keyword 'txt_float_default': PDF txt output floating-point format: default
struct txt_float_default_info {
  static const char* name() { return "default"; }
  static const char* help() { return "Use default formatting for "
  "floating-point output."; }
};
using txt_float_default = keyword< txt_float_default_info, d,e,f,a,u,l,t >;

// Keyword 'txt_float_scientific': PDF txt output floating-point format:
// scientific
struct txt_float_scientific_info {
  static const char* name() { return "scientific"; }
  static const char* help() { return "Use scientific formatting for "
  "floating-point output."; }
};
using txt_float_scientific =
  keyword< txt_float_scientific_info, s,c,i,e,n,t,i,f,i,c >;

// Keyword 'fixed': PDF txt output floating-point format: fixed
struct txt_float_fixed_info {
  static const char* name() { return "fixed"; }
  static const char* help() { return "Use fixed formatting for "
  "floating-point output."; }
};
using txt_float_fixed = keyword< txt_float_fixed_info, f,i,x,e,d >;

// Keyword 'precision': precision in digits for text output of numbers
struct precision_info {
  static const char* name() { return "precision"; }
  static const char* help() { return "Precision in digits for text output.."; }
};
using precision = keyword< precision_info, p,r,e,c,i,s,i,o,n >;

// Keyword 'centering': selected PDF output file centering
struct centering_info {
  static const char* name() { return "PDF output file centering"; }
  static const char* help() { return "PDF output file centering."; }
};
using pdf_centering = keyword< centering_info, c,e,n,t,e,r,i,n,g >;

// Keyword 'elem': PDF output file centering: elem
struct elem_info {
  static const char* name() { return "elem"; }
  static const char* help() { return "PDF output bins at elem centers."; }
};
using elem = keyword< elem_info, e,l,e,m >;

// Keyword 'node': PDF output file centering: node
struct node_info {
  static const char* name() { return "node"; }
  static const char* help() { return "PDF output bins at nodes."; }
};
using node = keyword< node_info, n,o,d,e >;


// Keyword 'init'
using init = keyword<undefined_info,  i,n,i,t >;

// Keyword 'coeff'
using coeff = keyword<undefined_info,  c,o,e,f,f >;

// Keyword 'raw': raw initialization policy
struct raw_info {
  static const char* name() { return "R"; }
  static const char* help() { return "Raw initialization policy."; }
};
using raw = keyword<raw_info, r,a,w >;

// Keyword 'raw': zero initialization policy
struct zero_info {
  static const char* name() { return "Z"; }
  static const char* help() { return "Zero initialization policy."; }
};
using zero = keyword<zero_info, z,e,r,o >;

// Keyword 'constant': constant coefficients policy
struct constant_info {
  static const char* name() { return "C"; }
  static const char* help() { return "Constant coefficients policy."; }
};
using constant = keyword<constant_info, c,o,n,s,t,a,n,t >;

// Keyword 'hommix'
struct hommix_info {
  static const char* name() { return "Homogeneous material mixing"; }
  static const char* help() { return
    "Physics option, 'hommix', is short for homogeneous material mixing. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test material mixing models independent, i.e., decoupled from other "
    "equations. Only a set of scalar equations are advanced which can be "
    "coupled to each other. The keyword 'hommix' introduces the hommix ... end "
    "block, selecting and describing the parameters of the mixing model(s). "
    "The common physics keywords are recognized.";
  }
};
using hommix = keyword<hommix_info, h,o,m,m,i,x>;

// Keyword 'homhydro'
struct homhydro_info {
  static const char* name() { return "Homogeneous hydrodynamics"; }
  static const char* help() { return
    "Physics option, 'homhydro', is short for homogeneous hydrodynamics. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test hydrodynamics models independent of, i.e., decoupled from other "
    "equations. Only a set of momentum equations are advanced whose components "
    "can be coupled to each other. The keyword 'homhydro' introduces the "
    "homhydro ... end block, selecting and describing the parameters of the "
    "hydrodynamics model(s). The common physics keywords are recognized.";
  }
};
using homhydro = keyword<homhydro_info,  h,o,m,h,y,d,r,o >;

// Keyword 'homrt'
struct homrt_info {
  static const char* name() { return "Homogeneous Rayleigh-Taylor"; }
  static const char* help() { return
    "Physics option, 'homrt', is short for homogeneous Rayleigh-Taylor. It is "
    "the simplest physics option that can be used to research, develop, and "
    "test hydrodynamics models for variable-density hydrodynamics and coupled "
    "material mixing, independent, i.e., decoupled from other equations. Only "
    "a set of mass and momentum conservation equations are advanced whose "
    "components can be coupled to each other. The keyword 'homrt' introduces "
    "the homrt ... end block, selecting and describing the parameters of the "
    "mass and hydrodynamics model(s). The common physics keywords are "
    "recognized.";
  }
};
using homrt = keyword<homrt_info,  h,o,m,r,t >;

// Keyword 'spinsflow'
struct spinsflow_info {
  static const char* name() { return
    "Standalone-particle incompressible Navier-Stokes flow";
  }
  static const char* help() { return
    "Physics option, 'spinsflow', is short for standalone-particle "
    "incompressible Navier-Stokes flow. It is a physics option intended for "
    "inhomogeneous constant-density flow. The transport equations solved are "
    "the momentum and optionally, energy, and a set of scalars. The "
    "divergence-constraint is enforced by solving a Poisson equation and "
    "projection scheme. The keyword 'spinsflow' introduces the spinsflow ... "
    "end block, selecting and describing the parameters of the above transport "
    "equations and their models. The common physics keywords are recognized.";
  }
};
using spinsflow = keyword<spinsflow_info,  s,p,i,n,s,f,l,o,w >;

// Keyword 'walker'
struct walker_info {
  static const char* name() {
    return "Walker. Differential equations testbed."; }
  static const char* help() { return
    "Test ordinary and stochastic differential equations.";
  }
};
using walker = keyword< walker_info, w,a,l,k,e,r >;

// Keyword 'dirichlet'
struct dirichlet_info {
  static const char* name() {
    return "Dirichlet"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "the Dirichlet distribution. For more details, see "
    "http://dx.doi.org/10.1155/2013/842981.";
  }
};
using dirichlet = keyword<dirichlet_info,  d,i,r,i,c,h,l,e,t >;

// Keyword 'generalized_dirichlet'
struct gendir_info {
  static const char* name() {
    return "Generalized Dirichlet"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "Lochner's generalized Dirichlet distribution. For more details, see "
    "http://dx.doi.org/10.1063/1.4822416.";
  }
};
using gendir =
  keyword< gendir_info, g,e,n,e,r,a,l,i,z,e,d,'_',d,i,r,i,c,h,l,e,t >;

// Keyword 'wright_fisher'
struct wrightfisher_info {
  static const char* name() {
    return "Wright-Fisher"; }
  static const char* help() { return
    "A coupled system of stochastic differential equations whose invariant is "
    "the Dirichlet distribution. For more details, see "
    "http://www.sciencedirect.com/science/article/pii/S0040580912001013.";
  }
};
using wrightfisher =
  keyword< wrightfisher_info, w,r,i,g,h,t,'_',f,i,s,h,e,r >;

// Keyword 'skewnormal'
struct skewnormal_info {
  static const char* name() {
    return "Skew-normal"; }
  static const char* help() { return
    "A system of stochastic differential equations whose invariant is the "
    " joint skew-normal distribution.";
  }
};
using skewnormal = keyword< skewnormal_info,  s,k,e,w,'-',n,o,r,m,a,l >;

// Keyword 'beta'
struct beta_info {
  static const char* name() {
    return "Beta"; }
  static const char* help() { return
    "A system of stochastic differential equations with linear drift and "
    "quadratic diagonal diffusion whose invariant is the joint beta "
    "distribution.";
  }
};
using beta = keyword< beta_info, b,e,t,a >;

// Keyword 'gamma'
struct gamma_info {
  static const char* name() {
    return "Gamma"; }
  static const char* help() { return
    "A system stochastic differential equations with linear drift and linear "
    "diagonal diffusion whose invariant is the joint gamma distribution.";
  }
};
using gamma = keyword< gamma_info, g,a,m,m,a >;

// Keyword 'ornstein_uhlenbeck'
struct ornstein_uhlenbeck_info {
  static const char* name() {
    return "Ornstein-Uhlenbeck"; }
  static const char* help() { return
    "A system of stochastic differential equations with linear drift and "
    "constant diffusion whose invariant is the joint normal distribution.";
  }
};
using ornstein_uhlenbeck =
  keyword< ornstein_uhlenbeck_info, o,r,n,s,t,e,i,n,'-',u,h,l,e,n,b,e,c,k >;

// Keyword 'diag_ornstein_uhlenbeck'
struct diag_ornstein_uhlenbeck_info {
  static const char* name() {
    return "Diagonal Ornstein-Uhlenbeck"; }
  static const char* help() { return
    "A system of stochastic differential equations with linear drift and "
    "constant diagonal diffusion whose invariant is the joint normal "
    "distribution.";
  }
};
using diag_ornstein_uhlenbeck =
  keyword< diag_ornstein_uhlenbeck_info,
           d,i,a,g,'_',o,r,n,s,t,e,i,n,'-',u,h,l,e,n,b,e,c,k >;

// Select position model:
//   * Insviscid model
using pos_inviscid = keyword<undefined_info,  p,o,s,'_',i,n,v,i,s,c,i,d >;
//   * Viscous model
using pos_viscous = keyword<undefined_info,  p,o,s,'_',v,i,s,c,o,u,s >;

// Select mass model:
//   * Beta model
using mass_beta = keyword<undefined_info,  m,a,s,s,'_',b,e,t,a >;

// Select hydrodynamics model:
//   * Simplified Langevin model
using hydro_slm = keyword<undefined_info,  h,y,d,r,o,'_',s,l,m >;
//   * Generalized Langevin model
using hydro_glm = keyword<undefined_info,  h,y,d,r,o,'_',g,l,m >;

// Keyword 'mix_iem'
struct mix_iem_info {
  static const char* name() { return "Interaction by exchange with the mean"; }
  static const char* help() { return
    "Material mix model, 'mix_iem', is short for interaction by exchange with "
    "the mean (IEM). It is a relaxation-type material mix model intended "
    "shear-driven flows.";
  }
};
using mix_iem = keyword<mix_iem_info,  m,i,x,'_',i,e,m >;

// Keyword 'mix_iecm'
struct mix_iecm_info {
  static const char* name() { return
    "Interaction by exchange with the conditional mean";
  }
  static const char* help() { return
    "Material mix model, 'mix_iecm', is short for interaction by exchange with "
    "the conditional mean (IECM). It is a relaxation-type material mix model "
    "intended shear-driven flows.";
  }
};
using mix_iecm = keyword<mix_iecm_info,  m,i,x,'_',i,e,c,m >;

// Keyword 'mix_dir'
struct mix_dir_info {
  static const char* name() { return "Dirichlet"; }
  static const char* help() { return
    "Material mix model, 'mix_dir', is short for Dirichlet. It is a material "
    "mix model that explicitly satisfies the unit-sum requirement for all "
    "statistical samples.";
  }
};
using mix_dir = keyword<mix_dir_info,  m,i,x,'_',d,i,r >;

// Keyword 'mix_gendir'
struct mix_gendir_info {
  static const char* name() { return "Generalized Dirichlet"; }
  static const char* help() { return
    "Material mix model, 'mix_gendir', is short for Lochner's generalized "
    "Dirichlet. It is a material mix model that explicitly satisfies the "
    "unit-sum requirement for all statistical samples.";
  }
};
using mix_gendir = keyword<mix_gendir_info,  m,i,x,'_',g,e,n,d,i,r >;

// Select material mix rate model:
//   * Gamma distribution model
using mixrate_gamma = keyword<undefined_info,  m,i,x,r,a,t,e,'_',g,a,m,m,a >;

// Select turbulence frequency model:
//   * Gamma distribution model
using freq_gamma = keyword<undefined_info,  f,r,e,q,'_',g,a,m,m,a >;

// Number of time steps to take
using nstep = keyword<undefined_info,  n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword<undefined_info,  t, e, r, m >;

// Size of time step
using dt = keyword<undefined_info,  d,t >;

// Start of position model specification block
using position = keyword<undefined_info,  p,o,s,i,t,i,o,n >;

// Start of hydrodynamics model specification block
using hydro = keyword<undefined_info,  h,y,d,r,o >;

// Start of material mix model specification block
using mix = keyword<undefined_info,  m,i,x >;

// Number of particle position components
using nposition = keyword<undefined_info,  n,p,o,s,i,t,i,o,n >;
// Number of particle density components
using ndensity = keyword<undefined_info,  n,d,e,n,s,i,t,y >;
// Number of particle velocity components
using nvelocity = keyword<undefined_info,  n,v,e,l,o,c,i,t,y >;
// Number of particle scalar components
using nscalar = keyword<undefined_info,  n,s,c,a,l,a,r >;
// Number of particle turbulence frequency components
using nfreq = keyword<undefined_info,  n,f,r,e,q >;

// Number of components
using ncomp = keyword<undefined_info,  n,c,o,m,p >;

// Wright-Fisher SDE parameters
using sde_omega = keyword<undefined_info, o,m,e,g,a >;

// Dirichlet and generalized Dirichlet parameters
using sde_b = keyword<undefined_info,  b >;
using sde_S = keyword<undefined_info,  S >;
using sde_kappa = keyword<undefined_info,  k,a,p,p,a >;
using sde_c = keyword<undefined_info,  c >;

// Ornstein-Uhlenbeck SDE parameters
using sde_sigma = keyword<undefined_info, s,i,g,m,a >;
using sde_theta = keyword<undefined_info, t,h,e,t,a >;
using sde_mu = keyword<undefined_info, m,u >;

// time scale
using sde_T = keyword<undefined_info, T >;

// lambda
using sde_lambda = keyword<undefined_info, l,a,m,b,d,a >;

// Langevin model parameters
using SLM_C0 = keyword<undefined_info,  C,'0' >;

// Gamma frequency model parameters
using freq_gamma_C1 = keyword<undefined_info,  C,'1' >;
using freq_gamma_C2 = keyword<undefined_info,  C,'2' >;
using freq_gamma_C3 = keyword<undefined_info,  C,'3' >;
using freq_gamma_C4 = keyword<undefined_info,  C,'4' >;

// Beta model parameters
using Beta_At = keyword<undefined_info,  A,t >;

// Quantities
using transported_scalar = keyword<undefined_info,  Y >;
using transported_scalar_fluctuation = keyword<undefined_info,  y >;

using velocity_x = keyword<undefined_info,  U >;
using velocity_fluctuation_x = keyword<undefined_info,  u >;
using velocity_y = keyword<undefined_info,  V >;
using velocity_fluctuation_y = keyword<undefined_info,  v >;
using velocity_z = keyword<undefined_info,  W >;
using velocity_fluctuation_z = keyword<undefined_info,  w >;

using pressure = keyword<undefined_info,  P >;
using pressure_fluctuation = keyword<undefined_info,  p >;
  
using density = keyword<undefined_info,  R >;
using density_fluctuation = keyword<undefined_info,  r >;
  
// Total number of particles
using npar = keyword<undefined_info,  n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword<undefined_info,  t,t,y,i >;

// Dump (restart file) output interval
using dmpi = keyword<undefined_info,  d,m,p,i >;

// Glob output interval
using glbi = keyword<undefined_info,  g,l,b,i >;

// Output interval
using interval = keyword<undefined_info,  i,n,t,e,r,v,a,l >;

// Statistics
using statistics = keyword<undefined_info,  s,t,a,t,i,s,t,i,c,s >;

// PDFs
using pdfs = keyword<undefined_info, p,d,f,s >;

// RNG block
using rngs = keyword<undefined_info,  r,n,g,s >;

// RNG
using rng = keyword<undefined_info,  r,n,g >;

// Keyword 'init'
using init = keyword<undefined_info,  i,n,i,t >;

// Keyword 'coeff'
using coeff = keyword<undefined_info,  c,o,e,f,f >;

// Keyword 'depvar': dependent variable in equations
struct depvar_info {
  static const char* name() { return "Dependent variable"; }
  static const char* help() { return "Dependent variable in differential "
                                     "equations."; }
};
using depvar = keyword< depvar_info, d,e,p,v,a,r >;

// Keyword 'control', cmdline '--control' with alias '-c'
struct control_info {
  static const char* name() { return "control"; }
  static const char* help() { return
    "This option is used to define the control file containing user input.";
  }
};
using control = cmdline_keyword< control_info, c, c,o,n,t,r,o,l >;


// Keyword 'smallcrush'
struct smallcrush_info {
  static const char* name() { return "SmallCrush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'SmallCrush'. SmallCrush is a "
    "battery of relatively small number, O(10), of tests, defined in TestU01, "
    "a library for the empirical testing of random number generators. For more "
    "info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using smallcrush = keyword< smallcrush_info, s,m,a,l,l,c,r,u,s,h >;

// Keyword 'crush'
struct crush_info {
  static const char* name() { return "Crush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'Crush'. Crush is a suite of "
    "stringent statistical tests, O(100), defined in TestU01, a library for "
    "the empirical testing of random number generators. For more info, see "
    "http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using crush = keyword< crush_info, c,r,u,s,h >;

// Keyword 'bigcrush'
struct bigcrush_info {
  static const char* name() { return "BigCrush"; }
  static const char* help() { return
    "This keyword is used to introduce the description of the random number "
    "generator test suite, i.e., battery, 'BigCrush'. BigCrush is a "
    "suite of very stringent statistical tests, O(100), defined in TestU01, a "
    "library for the empirical testing of random number generators. For more "
    "info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.";
  }
};
using bigcrush = keyword< bigcrush_info, b,i,g,c,r,u,s,h >;

// Keyword 'verbose', cmdline '--verbose' with alias '-v'
struct verbose_info {
  static const char* name() { return "verbose"; }
  static const char* help() { return
    "This option is used to pick verbose output.";
  }
};
using verbose = cmdline_keyword< verbose_info, v, v,e,r,b,o,s,e >;


// Keyword 'input', cmdline '--input' with alias '-i'
struct input_info {
  static const char* name() { return "input"; }
  static const char* help() { return
    "This option is used to define the input file.";
  }
};
using input = cmdline_keyword< input_info, i, i,n,p,u,t >;

// Keyword 'output', cmdline '--output' with alias '-o'
struct output_info {
  static const char* name() { return "output"; }
  static const char* help() { return
    "This option is used to define the output file.";
  }
};
using output = cmdline_keyword< output_info, o, o,u,t,p,u,t >;

// Keyword 'virtualization', cmdline '--virtualization' with alias '-u'
struct virtualization_info {
  static const char* name() { return "virtualization"; }
  static const char* help() { return
    "This option is used to define the degree of virtualization "
    "(over-decomposition). The virtualization parameter, is a real number "
    "between 0.0 and 1.0, inclusive, which controls the degree of "
    "virtualization or over-decomposition. Independent of the value of "
    "virtualization the work is approximately evenly distributed among the "
    "available processing elements. For zero virtualization (no "
    "over-decomposition), the work is simply decomposed into "
    "total_work/numPEs, which yields the smallest number of Charm++ chares and "
    "the largest chunks of work units. The other extreme is unity "
    "virtualization, which decomposes the total work into the smallest size "
    "work units possible, yielding the largest number of Charm++ chares. "
    "Obviously, the optimum will be between 0.0 and 1.0, depending on the "
    "problem.";
  }
};
using virtualization = cmdline_keyword< virtualization_info,
                                        u, v,i,r,t,u,a,l,i,z,a,t,i,o,n >;

// Keyword 'pdf', cmdline '--pdf' with alias '-p'
using pdf = cmdline_keyword<undefined_info, p, p,d,f >;

// Keyword 'stat', cmdline '--stat' with alias '-s'
using stat = cmdline_keyword<undefined_info, s, s,t,a,t >;

// Number of time steps to take
using nstep = keyword<undefined_info,  n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword<undefined_info,  t, e, r, m >;

// Size of time step
using dt = keyword<undefined_info,  d,t >;

// Number of components
using ncomp = keyword<undefined_info,  n,c,o,m,p >;

// Wright-Fisher SDE parameters
using sde_omega = keyword<undefined_info, o,m,e,g,a >;

// Dirichlet and generalized Dirichlet parameters
using sde_b = keyword<undefined_info,  b >;
using sde_S = keyword<undefined_info,  S >;
using sde_kappa = keyword<undefined_info,  k,a,p,p,a >;
using sde_c = keyword<undefined_info,  c >;

// Ornstein-Uhlenbeck SDE parameters
using sde_sigma = keyword<undefined_info, s,i,g,m,a >;
using sde_theta = keyword<undefined_info, t,h,e,t,a >;
using sde_mu = keyword<undefined_info, m,u >;

// time scale
using sde_T = keyword<undefined_info, T >;

// lambda
using sde_lambda = keyword<undefined_info, l,a,m,b,d,a >;

// Total number of particles
using npar = keyword<undefined_info,  n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword<undefined_info,  t,t,y,i >;

// Output interval
using interval = keyword<undefined_info,  i,n,t,e,r,v,a,l >;

// Statistics
using statistics = keyword<undefined_info,  s,t,a,t,i,s,t,i,c,s >;

// RNG block
using rngs = keyword<undefined_info,  r,n,g,s >;

// RNG
using rng = keyword<undefined_info,  r,n,g >;

// Keyword 'glob', cmdline '--glob' with alias '-g'
using glob = cmdline_keyword<undefined_info, g, g,l,o,b >;

// Select position model:
//   * Insviscid model
using pos_inviscid = keyword<undefined_info,  p,o,s,'_',i,n,v,i,s,c,i,d >;
//   * Viscous model
using pos_viscous = keyword<undefined_info,  p,o,s,'_',v,i,s,c,o,u,s >;

// Select mass model:
//   * Beta model
using mass_beta = keyword<undefined_info,  m,a,s,s,'_',b,e,t,a >;

// Select hydrodynamics model:
//   * Simplified Langevin model
using hydro_slm = keyword<undefined_info,  h,y,d,r,o,'_',s,l,m >;
//   * Generalized Langevin model
using hydro_glm = keyword<undefined_info,  h,y,d,r,o,'_',g,l,m >;

// Number of time steps to take
using nstep = keyword<undefined_info,  n,s,t,e,p >;

// Terminate time stepping at this value
using term = keyword<undefined_info,  t, e, r, m >;

// Size of time step
using dt = keyword<undefined_info,  d,t >;

// Start of position model specification block
using position = keyword<undefined_info,  p,o,s,i,t,i,o,n >;

// Start of hydrodynamics model specification block
using hydro = keyword<undefined_info,  h,y,d,r,o >;

// Start of material mix model specification block
using mix = keyword<undefined_info,  m,i,x >;

// Number of particle position components
using nposition = keyword<undefined_info,  n,p,o,s,i,t,i,o,n >;
// Number of particle density components
using ndensity = keyword<undefined_info,  n,d,e,n,s,i,t,y >;
// Number of particle velocity components
using nvelocity = keyword<undefined_info,  n,v,e,l,o,c,i,t,y >;
// Number of particle scalar components
using nscalar = keyword<undefined_info,  n,s,c,a,l,a,r >;
// Number of particle turbulence frequency components
using nfreq = keyword<undefined_info,  n,f,r,e,q >;

// Number of components
using ncomp = keyword<undefined_info,  n,c,o,m,p >;

// Wright-Fisher SDE parameters
using sde_omega = keyword<undefined_info, o,m,e,g,a >;

// Dirichlet and generalized Dirichlet parameters
using sde_b = keyword<undefined_info,  b >;
using sde_S = keyword<undefined_info,  S >;
using sde_kappa = keyword<undefined_info,  k,a,p,p,a >;
using sde_c = keyword<undefined_info,  c >;

// Ornstein-Uhlenbeck SDE parameters
using sde_sigma = keyword<undefined_info, s,i,g,m,a >;
using sde_theta = keyword<undefined_info, t,h,e,t,a >;
using sde_mu = keyword<undefined_info, m,u >;

// time scale
using sde_T = keyword<undefined_info, T >;

// lambda
using sde_lambda = keyword<undefined_info, l,a,m,b,d,a >;

// Langevin model parameters
using SLM_C0 = keyword<undefined_info,  C,'0' >;

// Gamma frequency model parameters
using freq_gamma_C1 = keyword<undefined_info,  C,'1' >;
using freq_gamma_C2 = keyword<undefined_info,  C,'2' >;
using freq_gamma_C3 = keyword<undefined_info,  C,'3' >;
using freq_gamma_C4 = keyword<undefined_info,  C,'4' >;

// Beta model parameters
using Beta_At = keyword<undefined_info,  A,t >;

// Quantities
using transported_scalar = keyword<undefined_info,  Y >;
using transported_scalar_fluctuation = keyword<undefined_info,  y >;

using velocity_x = keyword<undefined_info,  U >;
using velocity_fluctuation_x = keyword<undefined_info,  u >;
using velocity_y = keyword<undefined_info,  V >;
using velocity_fluctuation_y = keyword<undefined_info,  v >;
using velocity_z = keyword<undefined_info,  W >;
using velocity_fluctuation_z = keyword<undefined_info,  w >;

using pressure = keyword<undefined_info,  P >;
using pressure_fluctuation = keyword<undefined_info,  p >;
  
using density = keyword<undefined_info,  R >;
using density_fluctuation = keyword<undefined_info,  r >;
  
// Total number of particles
using npar = keyword<undefined_info,  n,p,a,r >;

// TTY (screen) output interval
using ttyi = keyword<undefined_info,  t,t,y,i >;

// Dump (restart file) output interval
using dmpi = keyword<undefined_info,  d,m,p,i >;

// Glob output interval
using glbi = keyword<undefined_info,  g,l,b,i >;

// Output interval
using interval = keyword<undefined_info,  i,n,t,e,r,v,a,l >;

// Statistics
using statistics = keyword<undefined_info,  s,t,a,t,i,s,t,i,c,s >;

// PDFs
using pdfs = keyword<undefined_info, p,d,f,s >;

// RNG block
using rngs = keyword<undefined_info,  r,n,g,s >;

// RNG
using rng = keyword<undefined_info,  r,n,g >;

// Keyword 'init'
using init = keyword<undefined_info,  i,n,i,t >;

// Keyword 'coeff'
using coeff = keyword<undefined_info,  c,o,e,f,f >;

} // kw::

#endif // Keywords_h
