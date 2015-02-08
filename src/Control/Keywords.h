//******************************************************************************
/*!
  \file      src/Control/Keywords.h
  \author    J. Bakosi
  \date      Sat 07 Feb 2015 08:45:09 PM MST
  \copyright 2012-2015, Jozsef Bakosi.
  \brief     Definition of all keywords
  \details   This file contains the definition of all keywords, including those
    of command-line argument parsers as well as input, i.e., control, file
    parsers. All keywords are shared among all parsers of all executables.

    All keywords are case-sensitive.

    The information contained in this file is used to build data structures for
    on-screen help on the command-line arguments and control file keywords,
    available via the --help, --helpctr, and --helpkw command-line arguments.

    A note on design: Defining structs that have static member functions
    returning a std::string is a way of storing C++-style strings at
    compile-time (which is not possible in a straightforward manner at this
    time). This could also be done with C-style const char* as well. The
    '*_info' structs store these strings, which then is used to specialize the
    _kw::keyword_ template, defined in Control/Keyword.h. Specializing the
    _keyword_ template also requires a specification of the precise string of
    characters that make up a keyword eventually matched by the parsers. Since
    these are all template arguments, the construction of keywords, their help,
    as well as all grammars, are entirely assembled at compile-time. Since the
    '*_info' struct member functions are static, they can be called without
    instantiating an object and thus available at compile-time.

    The definition of an '*_info' struct _requires_ at least the name, short,
    and long descriptions, defined by member functions _name(),
    _shortDescription()_, and _longDescription()_, respectively. Everything else
    is optional. However, when adding a new keyword it is highly recommended to
    define all of the _optional_ members if they make sense for the given
    keyword. If an expect value type is also given, that can be hooked up into
    where it is used.

    The most general definition of a keyword is as follows:

    \code{.cpp}
      // Keyword info definition
      struct keyword_info {

        // Required very short name, usually a single word or (for e.g.
        // policies) a single character. This can be the keyword itself, but
        // does not have to be. This field is used as an id of the option or
        // setting.
        static std::string name() { return "Name"; }

        // Required short keyword description
        static std::string shortDescription() { return "Short description"; }

        // Required detailed keyword description. This returns a string literal,
        // since this is usually multi-line and it is less work to maintain this
        // way.
        static std::string longDescription() { return
          R"(Longer, possibly multi-line description of the keyword. Example
          usage of the keyword is welcome here. Don't worry about formatting:
          when this field is printed, extra spaces will be removed and line
          breaks will be inserted.)";
        }

        // Optional keyword alias. See also kw::Alias and
        // <tpl_build>/pegtl/src/pegtl/include/pegtl/constants.hh for examples
        // of what can be passed as template arguments (basically single
        // characters). The one below defines the character 'c' as the alias.
        // Aliases are single character long. The idea of an alias is to have a
        // long as well as a short keyword for the same functionality.
        // Currently, this is only hooked up for command-line arguments and not
        // for control-file keywords, which is intentional. Command-line
        // arguments are a lot less than control file keywords and are more
        // frequently typed by the user. Thus command-line argument aliases are
        // user-friendly. There are many control file keywords and aliases only
        // cause confusion. Defining an alias for a command-line argument
        // enables the command-line parser to match on '--longer_keyword' as
        // well as on '-c'. Depending on whether the alias typedef is defined
        // for a keyword or not, the correct grammar is automatically generated
        // at compile-time, matching on both the longer keyword as well as on
        // the alias. Defining an alias for a control file keyword can be done
        // but has no effect in a control file parser.
        using alias = Alias< c >;

        // Optional expected data for the keyword - bundled to struct expect.
        // This struct is entirely optional within a keyword definition.
        // However, if it is defined, it must at least define the static member
        // function description() which returns the description of the type the
        // keyword expects. It may also optionally define the following fields:
        // 
        //    - type - defining the expected type
        //    - lower - lower bound of the expected value
        //    - upper - upper bound of the expected value
        //    - choices - valid choices for the expected value
        //
        struct expect {

          // If this struct is defined, required expected type description, max
          // 10 characters long
          static std::string description() { return "int"; }

          // Optional expected type
          using type = std::streamsize;

          // Optional expected value lower bound
          static const type lower = 1;

          // Optional expected value upper bound
          static const type upper = std::numeric_limits< tk::real >::digits10 + 1;

          // Optional expected valid choices description, here giving
          // information on the expected type and the valid bounds. Note that
          // this can be any string, but if they exist, it is a good idea give
          // at least some information on the bounds, as below, since the bounds
          // are NOT displayed in the help for a keyword. This decision keeps
          // the bounds specifications generic since they can be any type. As a
          // result, the help structures, defined in HelpFactory.h, are simpler
          // as they do not have to be parameterized by the type of the bounds,
          // which greatly simplifies that code.
          static std::string choices() {
            return "integer [" + std::to_string(lower) + "..." +
                   std::to_string(upper) + ']';
          }

        };

      };

      // Keyword definition, passing the above info struct as the first
      // template argument, and the rest of the template arguments are the
      // characters of the keyword to be matched by the parser. Remember: all
      // keywords are case-sensitive. This one contrived example below defines
      // keyword 'kw', matching the keyword 'KeYwOrD'.
      using kw = keyword< keyword_info, K,e,Y,w,O,r,D >;
    \endcode
  \see Control/Keyword.h
  \see Control/HelpFactory.h
*/
//******************************************************************************
#ifndef Keywords_h
#define Keywords_h

#include <Types.h>
#include <Keyword.h>

//! Keywords used by all input deck and command line parsers
namespace kw {

using namespace pegtl::ascii;

struct title_info {
  static std::string name() { return "title"; }
  static std::string shortDescription() { return "Set analysis title"; }
  static std::string longDescription() { return
    R"(The analysis title may be specified in the input file using the 'title'
    keyword. The 'title' keyword must be followed by a double-quoted string
    specifying the analysis title. Example: title "Example problem".
    Specifying a title is optional.)";
  }
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using title = keyword< title_info, t,i,t,l,e >;

struct end_info {
  static std::string name() { return "end"; }
  static std::string shortDescription() { return "End of an input block"; }
  static std::string longDescription() { return
    R"(The end of a block is given by the 'end' keyword in the input file.
    Example: "rngs ... end".)";
  }
};
using end = keyword< end_info, e,n,d >;

struct help_info {
  static std::string name() { return "help"; }
  static std::string shortDescription() { return
    R"(Display one-liner help on all command-line arguments)"; }
  static std::string longDescription() { return
    R"(Get a short one-liner help on all command-line arguments from an
    executable. It also triggers the help from the Charm++ runtime system and in
    addition to that of the executable, it also lists command-line arguments
    from Converse Machine, Tracing, Load Balancer, Record/Replay, and Charm++
    command-line parameters.)";
  }
  using alias = Alias< h >;
};
using help = keyword< help_info, h,e,l,p >;

struct helpctr_info {
  static std::string name() { return "helpctr"; }
  static std::string shortDescription() { return
    "Display one-liner help on all control file keywords"; }
  static std::string longDescription() { return
    R"(This keyword can be used to get a short one-liner help on all control
    file keywords from an executable.)";
  }
  using alias = Alias< f >;
};
using helpctr = keyword< helpctr_info, h,e,l,p,c,t,r >;

struct helpkw_info {
  static std::string name() { return "helpkw"; }
  static std::string shortDescription() { return
    "Display verbose help on a single keyword"; }
  static std::string longDescription() { return
    R"(This keyword can be used to get a verbose help on a single command-line
    argument or control-file keyword (i.e., help on keyword) from an
    executable.)";
  }
  using alias = Alias< H >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using helpkw = keyword< helpkw_info, h,e,l,p,k,w >;

struct seed_info {
  static std::string name() { return "seed"; }
  static std::string shortDescription() { return
    "Set random number generator seed"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a seed for a random number generator
    Example: rngmkl_mcg31 seed 1234 end)";
  }
  struct expect {
    using type = unsigned int;
    static std::string description() { return "uint"; }
  };
};
using seed = keyword< seed_info, s,e,e,d >;

struct mkl_mcg31_info {
  static std::string name() { return "MKL MCG311"; }
  static std::string shortDescription() { return
    "Select Intel MKL MCG31 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_MCG31', a 31-bit multiplicative
    congruential random number generator, provided by Intel's Math Kernel
    Library (MKL). For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_mcg31 = keyword< mkl_mcg31_info, m,k,l,'_',m,c,g,'3','1' >;

struct mkl_r250_info {
  static std::string name() { return "MKL R250"; }
  static std::string shortDescription() { return
    "Select Intel MKL R250 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_R250', a generalized feedback
    shift register random number generator, provided by Intel's Math Kernel
    Library (MKL). For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_r250 = keyword< mkl_r250_info, m,k,l,'_',r,'2','5','0' >;

struct mkl_mrg32k3a_info {
  static std::string name() { return "MKL MRG32K3A"; }
  static std::string shortDescription() { return
   "Select Intel MKL MRG32K3A RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_MRG32K3A', a combined multiple
    recursive random number generator with two components of order 3,
    provided by Intel's Math Kernel Library (MKL). For more info on MKL see
    https://software.intel.com/en-us/articles/intel-math-kernel-library-
    documentation.)";
  }
};
using mkl_mrg32k3a =
  keyword< mkl_mrg32k3a_info, m,k,l,'_',m,r,g,'3','2',k,'3',a >;

struct mkl_mcg59_info {
  static std::string name() { return "MKL MCG59"; }
  static std::string shortDescription() { return
    "Select Intel MKL MCG59 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_MCG59', a 59-bit multiplicative
    congruential random number generator, provided by Intel's Math Kernel
    Library (MKL). For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};

using mkl_mcg59 = keyword< mkl_mcg59_info, m,k,l,'_',m,c,g,'5','9' >;

struct mkl_wh_info {
  static std::string name() { return "MKL WH"; }
  static std::string shortDescription() { return
    "Select Intel MKL WH RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_WH', a set of 273 Wichmann-Hill
    combined multiplicative congruential random number generators, provided
    by Intel's Math Kernel Library (MKL). For more info on MKL see https://
    software.intel.com/en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_wh = keyword< mkl_wh_info, m,k,l,'_',w,h >;

struct mkl_mt19937_info {
  static std::string name() { return "MKL MT19937"; }
  static std::string shortDescription() { return
    "Select Intel MKL MT19937 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_MT19937', a Mersenne Twister
    pseudorandom number generator, provided by Intel's Math Kernel Library
    (MKL). For more info on MKL see https://software.intel.com/en-us/articles/
    intel-math-kernel-library-documentation.)";
  }
};
using mkl_mt19937 =
  keyword< mkl_mt19937_info, m,k,l,'_',m,t,'1','9','9','3','7' >;

struct mkl_mt2203_info {
  static std::string name() { return "MKL MT2203"; }
  static std::string shortDescription() { return
    "Select Intel MKL MT2203 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_MT2203', a set of 6024 Mersenne
    Twister pseudorandom number generators, available in Intel's Math Kernel
    Library (MKL). For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_mt2203 = keyword< mkl_mt2203_info, m,k,l,'_',m,t,'2','2','0','3' >;

struct mkl_sfmt19937_info {
  static std::string name() { return "MKL SFMT19937"; }
  static std::string shortDescription() { return
    "Select Intel MKL SFMT19937 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_SFMT19937', a SIMD-oriented Fast
    Mersenne Twister pseudorandom number generator, provided by Intel's Math
    Kernel Library (MKL). For more info on MKL see https://software.intel.com/
    en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_sfmt19937 =
  keyword< mkl_sfmt19937_info, m,k,l,'_',s,f,m,t,'1','9','9','3','7' >;

struct mkl_sobol_info {
  static std::string name() { return "MKL SOBOL"; }
  static std::string shortDescription() { return
    "Select Intel MKL SOBOL RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_SOBOL', a 32-bit Gray code-based
    random number generator, producing low-discrepancy sequences for
    dimensions 1 ≤ s ≤ 40 with available user-defined dimensions, provided
    by Intel's Math Kernel Library (MKL). For more info on MKL see https://
    software.intel.com/en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_sobol = keyword< mkl_sobol_info, m,k,l,'_',s,o,b,o,l >;

struct mkl_niederr_info {
  static std::string name() { return "MKL NIEDERR"; }
  static std::string shortDescription() { return
   "Select Intel MKL NIEDERR RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_NIEDERR', a 32-bit Gray
    code-based random number generator, producing low-discrepancy sequences
    for dimensions 1 ≤ s ≤ 318 with available user-defined dimensions,
    provided by Intel's Math Kernel Library (MKL). For more info on MKL see
    https://software.intel.com/en-us/articles/intel-math-kernel-library-
    documentation.)";
  }
};
using mkl_niederr = keyword< mkl_niederr_info, m,k,l,'_',n,i,e,d,e,r,r >;

struct mkl_iabstract_info {
  static std::string name() { return "MKL IABSTRACT"; }
  static std::string shortDescription() { return
    "Select Intel MKL IABSTRACT RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_IABSTRACT', an abstract random
    number generator for integer arrays, provided by Intel's Math Kernel
    Library (MKL). For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_iabstract =
  keyword< mkl_iabstract_info, m,k,l,'_',i,a,b,s,t,r,a,c,t >;

struct mkl_dabstract_info {
  static std::string name() { return "MKL DABSTRACT"; }
  static std::string shortDescription() { return
    "Select Intel MKL DABSTRACT RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_DABSTRACT', an abstract random
    number generator for double-precision floating-point arrays, provided by
    Intel's Math Kernel Library (MKL). For more info on MKL see https://
    software.intel.com/en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_dabstract =
  keyword< mkl_dabstract_info, m,k,l,'_',d,a,b,s,t,r,a,c,t >;

struct mkl_sabstract_info {
  static std::string name() { return "MKL SABSTRACT"; }
  static std::string shortDescription() { return
    "Select Intel MKL SABSTRACT RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_SABSTRACT', an abstract random
    number generator for single-precision floating-point arrays, provided by
    Intel's Math Kernel Library (MKL). For more info on MKL see https://
    software.intel.com/en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_sabstract =
  keyword< mkl_sabstract_info, m,k,l,'_',s,a,b,s,t,r,a,c,t >;

struct mkl_nondeterm_info {
  static std::string name() { return "MKL NONDETERM"; }
  static std::string shortDescription() { return
    "Select Intel MKL NONDETERM RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_NONDETERM', a non-deterministic
    random number generator, provided by Intel's Math Kernel Library (MKL).
    For more info on MKL see https://software.intel.com/en-us/articles/intel-
    math-kernel-library-documentation.)";
  }
};
using mkl_nondeterm =
  keyword< mkl_nondeterm_info, m,k,l,'_',n,o,n,d,e,t,e,r,m >;

struct standard_info {
  static std::string name() { return "standard"; }
  static std::string shortDescription() { return
    "Select the standard algorithm for uniform RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the standard method used to generate
    uniform random numbers using the Intel Math Kernel Library (MKL) random
    number generators. Valid options are 'standard' and 'accurate'. For more
    info on MKL see https://software.intel.com/en-us/articles/intel-math-
    kernel-library-documentation.)";
  }
};
using standard = keyword< standard_info, s,t,a,n,d,a,r,d >;

struct accurate_info {
  static std::string name() { return "accurate"; }
  static std::string shortDescription() { return
    "Select the accurate algorithm for uniform RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the accurate method used to generate
    uniform random numbers using the Intel Math Kernel Library (MKL) random
    number generators. Valid options are 'standard' and 'accurate'. For more
    info on MKL see https://software.intel.com/en-us/articles/intel-math-
    kernel-library-documentation.)";
  }
};
using accurate = keyword< accurate_info, a,c,c,u,r,a,t,e >;

struct uniform_method_info {
  static std::string name() { return "uniform method"; }
  static std::string shortDescription() { return
    "Select an Intel MKL uniform RNG method"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the method used to generate uniform
    random numbers using the Intel Math Kernel Library (MKL) random number
    generators. Valid options are 'standard' and 'accurate'. For more info on
    MKL see https://software.intel.com/en-us/articles/intel-math-kernel-
    library-documentation.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + standard::string() + "\' | \'"
                  + accurate::string() + '\'';
    }
  };
};
using uniform_method =
  keyword< uniform_method_info, u,n,i,f,o,r,m,'_',m,e,t,h,o,d >;

struct boxmuller_info {
  static std::string name() { return "Box-Muller"; }
  static std::string shortDescription() { return
   "Select the Box-Muller algorithm for sampling a Gaussian"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Box-Muller method used to generate
    Gaussian random numbers using the Intel Math Kernel Library (MKL) random
    random number generators. Valid options are 'boxmuller', 'boxmuller2',
    and 'icdf'. For more info on MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using boxmuller = keyword< boxmuller_info, b,o,x,m,u,l,l,e,r >;

struct boxmuller2_info {
  static std::string name() { return "Box-Muller 2"; }
  static std::string shortDescription() { return
   "Select the Box-Muller 2 algorithm for sampling a Gaussian"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the Box-Muller 2 method used to generate
    Gaussian random numbers using the Intel Math Kernel Library (MKL) random
    number generators. For more info on MKL see https://software.intel.com/
    en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using boxmuller2 = keyword< boxmuller2_info, b,o,x,m,u,l,l,e,r,'2' >;

struct icdf_info {
  static std::string name() { return "ICDF"; }
  static std::string shortDescription() { return
    R"(Use inverse cumulative distribution function for sampling a Gaussian)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the inverse cumulative distribution
    function (ICDF) method used to generate Gaussian random numbers using the
    Intel Math Kernel Library (MKL) random number generators. For more info
    on MKL see https://software.intel.com/en-us/articles/intel-math-kernel-
    library-documentation.)";
  }
};
using icdf = keyword< icdf_info, i,c,d,f >;

struct gaussian_method_info {
  static std::string name() { return "Gaussian method"; }
  static std::string shortDescription() { return
    "Select an Intel MKL Gaussian RNG method"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the method used to generate Gaussian
    random numbers using the Intel Math Kernel Library (MKL) random number
    generators. Valid options are 'boxmuller', 'boxmuller2', and 'icdf'. For
    more info on MKL see https://software.intel.com/en-us/articles/intel-math-
    kernel-library-documentation.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + boxmuller::string() + "\' | \'"
                  + boxmuller2::string() + "\' | \'"
                  + icdf::string() + '\'';
    }
  };
};
using gaussian_method =
  keyword< gaussian_method_info, g,a,u,s,s,i,a,n,'_',m,e,t,h,o,d >;

struct rngsse_gm19_info {
  static std::string name() { return "RNGSSE GM19"; }
  static std::string shortDescription() { return
    "Select RNGSSE GM19 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GM19 random number generator, using a
    method based on parallel evolution of an ensemble of transformations of
    a two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gm19 = keyword< rngsse_gm19_info, r,n,g,s,s,e,'_',g,m,'1','9' >;

struct rngsse_gm29_info {
  static std::string name() { return "RNGSSE GM29"; }
  static std::string shortDescription() { return
    "Select RNGSSE GM29 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GM29 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gm29 = keyword< rngsse_gm29_info, r,n,g,s,s,e,'_',g,m,'2','9' >;

struct rngsse_gm31_info {
  static std::string name() { return "RNGSSE GM31"; }
  static std::string shortDescription() { return
    "Select RNGSSE GM31 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GM31 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gm31 = keyword< rngsse_gm31_info, r,n,g,s,s,e,'_',g,m,'3','1' >;

struct rngsse_gm55_info {
  static std::string name() { return "RNGSSE GM55"; }
  static std::string shortDescription() { return
    "Select RNGSSE GM55 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GM55 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gm55 = keyword< rngsse_gm55_info, r,n,g,s,s,e,'_',g,m,'5','5' >;

struct rngsse_gm61_info {
  static std::string name() { return "RNGSSE GM61"; }
  static std::string shortDescription() { return
    "Select RNGSSE GM61 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GM61 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gm61 = keyword< rngsse_gm61_info, r,n,g,s,s,e,'_',g,m,'6','1' >;

struct rngsse_gq581_info {
  static std::string name() { return "RNGSSE GQ58.1"; }
  static std::string shortDescription() { return
    "Select RNGSSE GQ58.1 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GQ58.1 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gq581 =
  keyword< rngsse_gq581_info, r,n,g,s,s,e,'_',g,q,'5','8','.','1' >;

struct rngsse_gq583_info {
  static std::string name() { return "RNGSSE GQ58.3"; }
  static std::string shortDescription() { return
    "Select RNGSSE GQ58.3 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GQ58.3 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSS2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gq583 =
  keyword< rngsse_gq583_info, r,n,g,s,s,e,'_',g,q,'5','8','.','3' >;

struct rngsse_gq584_info {
  static std::string name() { return "RNGSSE GQ58.4"; }
  static std::string shortDescription() { return
    "Select RNGSSE GQ58.4 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the GQ58.4 generator, using a method based
    on parallel evolution of an ensemble of transformations of a
    two-dimensional torus, provided by the RNGSSE2 random number generator
    library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_gq584 =
  keyword< rngsse_gq584_info, r,n,g,s,s,e,'_',g,q,'5','8','.','4' >;

struct rngsse_mt19937_info {
  static std::string name() { return "RNGSSE MT19937"; }
  static std::string shortDescription() { return
    "Select RNGSSE MT19937 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the MT19937 generator, a Mersenne Twister
    generator, provided by the RNGSSE2 random number generator library. For
    more info on RNGSSE see http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_mt19937 =
  keyword< rngsse_mt19937_info, r,n,g,s,s,e,'_',m,t,'1','9','9','3','7' >;

struct rngsse_lfsr113_info {
  static std::string name() { return "RNGSSE LSFR113"; }
  static std::string shortDescription() { return
    "Select RNGSSE LFSR113 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the LFSR113 generator, provided by the
    RNGSSE2 random number generator library. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_lfsr113 =
  keyword< rngsse_lfsr113_info, r,n,g,s,s,e,'_',l,f,s,r,'1','1','3' >;

struct rngsse_mrg32k3a_info {
  static std::string name() { return "RNGSSE MRG32K3A"; }
  static std::string shortDescription() { return
    "Select RNGSSE MRG32K3A RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the MRG32K3A generator, a combined
    multiple recursive random number generator with two components of order
    3, provided by the RNGSS2 random number generator library. For more info
    on RNGSSE see http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using rngsse_mrg32k3a =
  keyword< rngsse_mrg32k3a_info, r,n,g,s,s,e,'_',m,r,g,'3','2',k,'3',a >;

struct seq_short_info {
  static std::string name() { return "short"; }
  static std::string shortDescription() { return
    "Select the short sequence length for an RNGSSE2 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the short sequence length used by the
    RNGSSE2 random number generator library. Valid options are 'short',
    'medium', and 'long'. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using seq_short = keyword< seq_short_info, s,h,o,r,t >;

struct seq_med_info {
  static std::string name() { return "medium"; }
  static std::string shortDescription() { return
    "Select the medium sequence length for an RNGSSE2 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the medium sequence length used by the
    RNGSSE2 random number generator library. Valid options are 'short',
    'medium', and 'long'. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using seq_med = keyword< seq_med_info, m,e,d,i,u,m >;

struct seq_long_info {
  static std::string name() { return "long"; }
  static std::string shortDescription() { return
    "Select the long sequence length for an RNGSSE2 RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the medium sequence length used by the
    RNGSSE2 random number generator library. Valid options are 'short',
    'medium', and 'long'. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
};
using seq_long = keyword< seq_long_info, l,o,n,g >;

struct seqlen_info {
  static std::string name() { return "RNGSSE2 sequence length"; }
  static std::string shortDescription() { return
    "Specify the RNGSSE RNG sequence length"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a random number generator sequence length,
    used by the RNGSSE2 random number generator library. Valid options are
    'short', 'medium', and 'long'. For more info on RNGSSE see
    http://dx.doi.org/10.1016/j.cpc.2011.03.022.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + seq_short::string() + "\' | \'"
                  + seq_med::string() + "\' | \'"
                  + seq_long::string() + '\'';
    }
  };
};
using seqlen = keyword< seqlen_info, s,e,q,l,e,n >;

struct pdfs_info {
  static std::string name() { return "PDFs block"; }
  static std::string shortDescription() { return
    "Start of probability density function (PDF) input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested output for probability density
    functions (PDFs). Example: "pdfs mypdf( y1 : 1.0e-2 ) end", which
    requests a single-variate PDF to be output to file, whose sample space
    variable is y1, using automatic determination of the bounds of the sample
    space, using 1.0e-2 as the sample space bin size, and call the PDF
    "mypdf". For more info on the structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox.)";
  }
};
using pdfs = keyword< pdfs_info, p,d,f,s >;

struct txt_info {
  static std::string name() { return "txt"; }
  static std::string shortDescription() { return
    "Select ASCII output for outputing PDFs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the text
    output file type of a requested probability density function (PDF) within
    a pdfs ... end block. Example: "filetype txt", which selects text-file
    output. Valid options are 'txt', 'gmshtxt', 'gmshbin', and 'exodusii'.
    For more info on the structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox.)"; }
};
using txt = keyword< txt_info, t,x,t >;

struct gmshtxt_info {
  static std::string name() { return "gmshtxt"; }
  static std::string shortDescription() { return
    "Select Gmsh ASCII output for outputing PDFs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the ASCII
    (text) output file type readable by Gmsh of a requested probability
    density function (PDF) within a pdfs ... end block. Example: "filetype
    gmshtxt", which selects Gmsh ASCII file output. Valid options are 'txt',
    'gmshtxt', 'gmshbin', and 'exodusii'. For more info on the structure of
    the pdfs ... end block, see doc/pages/statistics_output.dox. For more
    info on Gmsh, see http://www.geuz.org/gmsh.)"; }
};
using gmshtxt = keyword< gmshtxt_info, g,m,s,h,t,x,t >;

struct gmshbin_info {
  static std::string name() { return "gmshbin"; }
  static std::string shortDescription() { return
    "Select Gmsh binary output for outputing PDFs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    binary output file type readable by Gmsh of a requested probability
    density function (PDF) within a pdfs ... end block. Example: "filetype
    gmshbin", which selects Gmsh binary file output. Valid options are 'txt',
    'gmshtxt', 'gmshbin', and 'exodusii'. For more info on the structure of
    the pdfs ... end block, see doc/pages/statistics_output.dox. For more
    info on Gmsh, see http://www.geuz.org/gmsh.)"; }
};
using gmshbin = keyword< gmshbin_info, g,m,s,h,b,i,n >;

struct exodusii_info {
  static std::string name() { return "exo"; }
  static std::string shortDescription() { return
    "Select ExodusII output for outputing PDFs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    ExodusII output file type readable by, e.g., ParaView of a requested
    probability density function (PDF) within a pdfs ... end block. Example:
    "filetype exodusii", which selects ExodusII file output. Valid options
    are 'txt', 'gmshtxt', 'gmshbin', and 'exodusii'. For more info on the
    structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox. For more info on ExodusII, see
    http://sourceforge.net/projects/exodusii.)";
  }
};
using exodusii = keyword< exodusii_info, e,x,o,d,u,s,i,i >;

struct filetype_info {
  static std::string name() { return "filetype"; }
  static std::string shortDescription() { return
    "Select PDF output file type"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the output file type of a requested
    probability density function (PDF) within a pdfs ... end block. Example:
    "filetype exodusii", which selects ExodusII output. Valid options are
    'txt', 'gmshtxt', 'gmshbin', and 'exodusii'. For more info on the
    structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + txt::string() + "\' | \'"
                  + gmshtxt::string() + "\' | \'"
                  + gmshbin::string() + "\' | \'"
                  + exodusii::string() + '\'';
    }
  };

};
using pdf_filetype = keyword< filetype_info, f,i,l,e,t,y,p,e >;

struct overwrite_info {
  static std::string name() { return "overwrite"; }
  static std::string shortDescription() { return
    "Select PDF output policy overwrite"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    the 'overwrite' output file policy for requested probability density
    functions (PDFs) within a pdfs ... end block. Example: "policy
    overwrite", which selects the overwrite output file policy. The
    overwrite policy overwrites the same output file containing a single time
    step. Valid PDF policy options are 'overwrite', 'multiple', and
    'evolution'. For more info on the structure of the pdfs ... end block,
    see doc/pages/statistics_output.dox.)"; }
};
using overwrite = keyword< overwrite_info, o,v,e,r,w,r,i,t,e >;

struct multiple_info {
  static std::string name() { return "multiple"; }
  static std::string shortDescription() { return
    "Select PDF output policy multiple"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    the 'multiple' output file policy for requested probability density
    functions (PDFs) within a pdfs ... end block. Example: "policy
    multiple", which selects the multiple output file policy. The
    multiple policy output creates a new file for each time step. Valid PDF
    policy options are 'overwrite', 'multiple', and 'evolution'. For more
    info on the structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox.)"; }
};
using multiple = keyword< multiple_info, m,u,l,t,i,p,l,e >;

struct evolution_info {
  static std::string name() { return "evolution"; }
  static std::string shortDescription() { return
    "Select PDF output policy evolution"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    the 'evolution' output file policy for requested probability density
    functions (PDFs) within a pdfs ... end block. Example: "policy
    evolution", which selects the evolution output file policy. The
    evolution policy output appends new time step to the same output file for
    each time instant, yielding a time evolution of data in a single file.
    Valid PDF policy options are 'overwrite', 'multiple', and 'evolution'. For
    more info on the structure of the pdfs ... end block, see
    doc/pages/statistics_output.dox.)";
  }
};
using evolution = keyword< evolution_info, e,v,o,l,u,t,i,o,n >;

struct policy_info {
  static std::string name() { return "policy"; }
  static std::string shortDescription() { return
    "Select PDF output file policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    the output file policy for requested probability density functions
    (PDFs) within a pdfs ... end block. Example: "policy overwrite", which
    selects the overwrite output file policy. Valid options are 'overwrite',
    'multiple', and 'evolution'. For more info on the structure of the
    pdfs ... end block, see doc/pages/statistics_output.dox.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + overwrite::string() + "\' | \'"
                  + multiple::string() + "\' | \'"
                  + evolution::string() + '\'';
    }
  };
};
using pdf_policy = keyword< policy_info, p,o,l,i,c,y >;

struct txt_float_default_info {
  static std::string name() { return "default"; }
  static std::string shortDescription() { return
   "Select the default ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'default' floating-point output format for ASCII floating-point real
    number output. Example: "format default", which selects the default
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_default = keyword< txt_float_default_info, d,e,f,a,u,l,t >;

struct txt_float_scientific_info {
  static std::string name() { return "scientific"; }
  static std::string shortDescription() { return
   "Select the scientific ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'scientific' floating-point output format for ASCII floating-point real
    number output. Example: "format scientific", which selects the
    scientific floating-point output. Valid options are 'default', 'fixed',
    and 'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_scientific =
  keyword< txt_float_scientific_info, s,c,i,e,n,t,i,f,i,c >;

struct txt_float_fixed_info {
  static std::string name() { return "fixed"; }
  static std::string shortDescription() { return
   "Select the fixed ASCII floating-point output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    'fixed' floating-point output format for ASCII floating-point real
    number output. Example: \"format fixed\", which selects the fixed
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)";
  }
};
using txt_float_fixed = keyword< txt_float_fixed_info, f,i,x,e,d >;

struct txt_float_format_info {
  static std::string name() { return "float format"; }
  static std::string shortDescription() { return
    "Specify the ASCII floating-point output format"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    floating-point output format for ASCII floating-point number output.
    Example: "format scientific", which selects the scientific
    floating-point output. Valid options are 'default', 'fixed', and
    'scientific'. For more info on these various formats, see
    http://en.cppreference.com/w/cpp/io/manip/fixed.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + txt_float_default::string() + "\' | \'"
                  + txt_float_scientific::string() + "\' | \'"
                  + txt_float_fixed::string() + '\'';
    }
  };
};
using txt_float_format = keyword< txt_float_format_info, f,o,r,m,a,t >;

struct precision_info {
  static std::string name() { return "precision"; }
  static std::string shortDescription() { return
    R"(Precision in digits for ASCII floating-point output)"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    the precision in digits for ASCII floating-point real number output.
    Example: "precision 10", which selects ten digits for floating-point
    output, e.g., 3.141592654. The number of digits must be larger than zero
    and lower than the maximum representable digits for the given floating-point
    type, defined by std::numeric_limits< FLOAT_TYPE >::digits10+2.
    For more info on setting the precision in C++, see
    http://en.cppreference.com/w/cpp/io/manip/setprecision, and
    http://en.cppreference.com/w/cpp/types/numeric_limits/digits10)";
  }
  struct expect {
    using type = std::streamsize;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< tk::real >::digits10 + 1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using precision = keyword< precision_info, p,r,e,c,i,s,i,o,n >;

struct elem_info {
  static std::string name() { return "elem"; }
  static std::string shortDescription() { return
    "Specify elem-centering for PDF output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    element-centering for the probability values on the sample space
    grid for file output of probability density functions (PDFs). Example:
    "centering elem", which selects element-centered values. Valid options
    are 'elem' and 'node', denoting cell-centered and point-centered output,
    respectively.)"; }
};
using elem = keyword< elem_info, e,l,e,m >;

struct node_info {
  static std::string name() { return "node"; }
  static std::string shortDescription() { return
    "Specify node-centering for PDF output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select
    node-centering for the probability values on the sample space grid for
    file output of probability density functions (PDFs). Example: "centering
    elem", which selects element-centered values. Valid options are 'elem'
    and 'node', denoting cell-centered and point-centered output,
    respectively.)"; }
};
using node = keyword< node_info, n,o,d,e >;

struct centering_info {
  static std::string name() { return "centering"; }
  static std::string shortDescription() { return
    "Specify data-centering for PDF output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the data
    centering of the probability value output on the sample space grid for
    file output of probability density functions (PDFs). Example: "centering
    elem", which selects element-centered values. Valid options are 'elem'
    and 'node', denoting cell-centered and point-centered output,
    respectively.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + elem::string() + "\' | \'"
                  + node::string() + '\'';
    }
  };
};
using pdf_centering = keyword< centering_info, c,e,n,t,e,r,i,n,g >;

struct raw_info {
  static std::string name() { return "R"; }
  static std::string shortDescription() { return
    "Select the raw initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the raw
    initialization policy. The initialization policy is used to specify how
    the initial conditions are set at t = 0 before time-integration.
    Example: "init raw", which selects raw initialization policy, which
    leaves the memory uninitialized. Note that this option may behave
    differently depending on the particular equation or physical model. For
    an example, see tk::InitPolicies in DiffEq/InitPolicy.h for valid
    options.)"; }
};
using raw = keyword< raw_info, r,a,w >;

struct zero_info {
  static std::string name() { return "Z"; }
  static std::string shortDescription() { return
    "Select the zero initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the zero
    initialization policy. The initialization policy is used to specify how
    the initial conditions are set at t = 0 before time-integration.
    Example: "init zero", which selects zero initialization policy, which
    puts zeros in memory. Note that this option may behave differently
    depending on the particular equation or physical model. For an example,
    see tk::InitPolicies in DiffEq/InitPolicy.h for valid options.)"; }
};
using zero = keyword< zero_info, z,e,r,o >;

struct init_info {
  static std::string name() { return "initpolicy"; }
  static std::string shortDescription() { return
    "Select initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select an
    initialization policy. This is used to specify how the initial conditions
    are set at t = 0 before time-integration. Example: "init raw", which
    selects raw initialization policy, which leaves the memory uninitialized.
    Note that this option may behave differently depending on the particular
    equation or physical model. For an example, see tk::InitPolicies in
    DiffEq/InitPolicy.h for valid options.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + raw::string() + "\' | \'"
                  + zero::string() + '\'';
    }
  };
};
using init = keyword< init_info,  i,n,i,t >;

struct const_info {
  static std::string name() { return "C"; }
  static std::string shortDescription() { return
    "Select constant coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    constant coefficients policy. The coefficients policy is used to specify
    how the coefficients are set at each time step during time-integration.
    Example: "coeff const", which selects constant coefficients policy,
    which sets constant coefficients before t = 0 and leaves the coefficients
    unchanged during time integration. Note that this option may behave
    differently depending on the particular equation or physical model. For
    an example, see walker::DirCoeffPolicies in DiffEq/DirCoeffPolicy.h for
    valid options.)"; }
};
using constant = keyword< const_info, c,o,n,s,t >;

struct jrrj_info {
  static std::string name() { return "J"; }
  static std::string shortDescription() { return
    "Select JRRJ coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the JRRJ coefficients policy. This policy
    (or model) is named after Joseph Raymond Ristorcelli Jr., which, at this
    time, is used to test some of Ray's closure ideas while we are working on
    a PDF/moment closure for variable-density binary material mixing for
    turbulent flows based on the beta distribution. A coefficients policy, in
    general, is used to specify how
    the coefficients are set at each time step during time-integration. Example:
    "coeff const", which selects constant coefficients policy, which sets
    constant coefficients before t = 0 and leaves the coefficients unchanged
    during time integration. Note that this option may behave differently
    depending on the particular equation or physical model. For an example, see
    walker::DirCoeffPolicies in DiffEq/DirCoeffPolicy.h for valid options.)"; }
};
using jrrj = keyword< jrrj_info, j,r,r,j >;

struct coeff_info {
  static std::string name() { return "coeffpolicy"; }
  static std::string shortDescription() { return
    "Select the coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a
    coefficients policy. This is used to specify how the coefficients are set
    at each time step during time-integration. Example: "coeff const",
    which selects constant coefficients policy, which sets constant
    coefficients before t = 0 and leaves the coefficients unchanged during
    time integration. Note that this option may behave differently depending
    on the particular equation or physical model. For an example, see
    walker::DirCoeffPolicies in DiffEq/DirCoeffPolicy.h for valid options.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + constant::string() + '\'';
    }
  };
};
using coeff = keyword< coeff_info,  c,o,e,f,f >;

struct walker_info {
  static std::string name() { return "walker"; }
  static std::string shortDescription() { return
    "Start configuration block of the random walker"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the walker. Walker, is a random walker,
    that allows temporal integration of a system of ordinary or stochastic
    differential equations (SDEs) of various types and the collection of
    arbitrary coupled statistics and probability density functions. Walker is
    intended as a general mathematical tool to analyze the behavior of SDEs
    and its statistics.)";
  }
};
using walker = keyword< walker_info, w,a,l,k,e,r >;

struct npar_info {
  static std::string name() { return "npar"; }
  static std::string shortDescription() { return
    "Set total number of particles"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the total number of particles in a
    simulation.)";
  }
  struct expect {
    using type = uint64_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using npar = keyword< npar_info, n,p,a,r >;

struct nstep_info {
  static std::string name() { return "nstep"; }
  static std::string shortDescription() { return
    "Set number of time steps to take"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of time steps to take in a
    simulation. The number of time steps are used in conjunction with the
    maximmum time specified by keyword 'term': the simulation stops whichever is
    reached first. Both 'nstep' and 'term' can be left unspecified, in which
    case their default values are used. See also 'term'.)";
  }
  struct expect {
    using type = uint64_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using nstep = keyword< nstep_info, n,s,t,e,p >;

struct term_info {
  static std::string name() { return "term"; }
  static std::string shortDescription() { return
    "Set maximum non-dimensional time to simulate"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the termination time in a simulation. The
    termination time and number of time steps, specified by 'nstep', are used in
    conjunction to determine when to stop a simulation: whichever is
    reached first. Both 'nstep' and 'term' can be left unspecified, in which
    case their default values are used. See also 'nstep'.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using term = keyword< term_info, t, e, r, m >;

struct dt_info {
  static std::string name() { return "dt"; }
  static std::string shortDescription() { return
    "Select (initial) time step size"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the time step size. For a variable-dt
    simulation this is the initial time step size.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using dt = keyword< dt_info, d,t >;

struct ncomp_info {
  static std::string name() { return "ncomp"; }
  static std::string shortDescription() { return
    "Set number of scalar components for a system of SDEs"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of scalar components of a
    vector.)";
  }
  struct expect {
    using type = unsigned int;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using ncomp = keyword< ncomp_info,  n,c,o,m,p >;

struct ttyi_info {
  static std::string name() { return "ttyi"; }
  static std::string shortDescription() { return
    "Set screen output interval"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the interval in time steps for screen
    output during a simulation.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using ttyi = keyword< ttyi_info, t,t,y,i >;

struct interval_info {
  static std::string name() { return "interval"; }
  static std::string shortDescription() { return
    "Set interval (within a relevant block)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify an interval in time steps. This must be
    used within a relevant block.)";
  }
  struct expect {
    using type = uint32_t;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using interval = keyword< interval_info, i,n,t,e,r,v,a,l >;

struct statistics_info {
  static std::string name() { return "statistics"; }
  static std::string shortDescription() { return
    "Start of statistics input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested output for statistical moments.
    Example: "statistics <Y> <yy> end", which requests the first two moments of
    the flutcutating variable 'Y'. For more info on the structure of the
    statistics ... end block, see doc/pages/statistics_output.dox.)";
  }
};
using statistics = keyword< statistics_info, s,t,a,t,i,s,t,i,c,s >;

struct rngs_info {
  static std::string name() { return "rngs"; }
  static std::string shortDescription() { return
    "Start of a random number generators description input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    descriptions and settings of requested random number generators.
    Example: "rngs mkl_mcg59 seed 2134 uniform_method accurate end end" which
    enables the MCG59 generator from MKL using the seed 2134. For more info on
    the structure of the rngs ... end block, see
    doc/pages/rngs_input.dox.)";
  }
};
using rngs = keyword< rngs_info, r,n,g,s >;

struct rng_info {
  static std::string name() { return "rng"; }
  static std::string shortDescription() { return
    "Select random number generator (RNG) from pool of enabled RNGs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a particular random number generator (RNG)
    from a pre-selected set of (enabled and configured) pool of RNGs. The pool
    is specified by the 'rngs ... end' block and it must precede the selection
    of an RNG.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + rngsse_gm19::string()     + "\' | \'"
                  + rngsse_gm29::string()     + "\' | \'"
                  + rngsse_gm31::string()     + "\' | \'"
                  + rngsse_gm55::string()     + "\' | \'"
                  + rngsse_gm61::string()     + "\' | \'"
                  + rngsse_gq581::string()    + "\' | \'"
                  + rngsse_gq583::string()    + "\' | \'"
                  + rngsse_gq584::string()    + "\' | \'"
                  + rngsse_mt19937::string()  + "\' | \'"
                  + rngsse_lfsr113::string()  + "\' | \'"
                  + rngsse_mrg32k3a::string()
                  #ifdef HAS_MKL
                                              + "\' | \'"
                  + mkl_mcg31::string()       + "\' | \'"
                  + mkl_r250::string()        + "\' | \'"
                  + mkl_mrg32k3a::string()    + "\' | \'"
                  + mkl_mcg59::string()       + "\' | \'"
                  + mkl_wh::string()          + "\' | \'"
                  + mkl_mt19937::string()     + "\' | \'"
                  + mkl_mt2203::string()      + "\' | \'"
                  + mkl_sfmt19937::string()   + "\' | \'"
                  + mkl_sobol::string()       + "\' | \'"
                  + mkl_niederr::string()     + "\' | \'"
                  + mkl_iabstract::string()   + "\' | \'"
                  + mkl_dabstract::string()   + "\' | \'"
                  + mkl_sabstract::string()   + "\' | \'"
                  + mkl_nondeterm::string()   + "\' "
                  #else
                                              + "\' "
                  #endif
             + "Remember: the RNG must be listed in the pool before it can be "
               "selected via this keyword!";
    }
  };
};
using rng = keyword< rng_info, r,n,g >;

struct sde_omega_info {
  static std::string name() { return "omega"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) omega)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "omega 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_omega = keyword< sde_omega_info, o,m,e,g,a >;

struct sde_b_info {
  static std::string name() { return "b"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) b)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "b 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_b = keyword< sde_b_info,  b >;

struct sde_S_info {
  static std::string name() { return "S"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) S)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "S 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_S = keyword< sde_S_info,  S >;

struct sde_kappa_info {
  static std::string name() { return "kappa"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) kappa)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "kappa 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_kappa = keyword< sde_kappa_info,  k,a,p,p,a >;

struct sde_c_info {
  static std::string name() { return "c"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) c)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "c 5.0 2.0 3.0 end". The length of the vector depends on the particular type
    of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_c = keyword< sde_c_info,  c >;

struct sde_sigmasq_info {
  static std::string name() { return "sigmasq"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) sigmasq)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "sigmasq 5.0 2.0 3.0 end". The length of the vector depends on the
    particular type of SDE system and is controlled by the preceding keyword
    'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_sigmasq = keyword< sde_sigmasq_info, s,i,g,m,a,s,q >;

struct sde_theta_info {
  static std::string name() { return "theta"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) theta)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "theta 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_theta = keyword< sde_theta_info, t,h,e,t,a >;

struct sde_mu_info {
  static std::string name() { return "mu"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) mu)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "mu 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_mu = keyword< sde_mu_info, m,u >;

struct sde_T_info {
  static std::string name() { return "T"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) T)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "T 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_T = keyword< sde_T_info, T >;

struct sde_lambda_info {
  static std::string name() { return "lambda"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) lambda)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "lambda 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_lambda = keyword< sde_lambda_info, l,a,m,b,d,a >;

struct depvar_info {
  static std::string name() { return "depvar"; }
  static std::string shortDescription() { return
    "Select dependent variable (in a relevant block)"; }
  static std::string longDescription() { return
    R"(Dependent variable, e.g, in differential equations.)"; }
};
using depvar = keyword< depvar_info, d,e,p,v,a,r >;

struct dirichlet_info {
  static std::string name() { return "Dirichlet"; }
  static std::string shortDescription() { return
    "Start configuration block for the Dirichlet SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the dirichlet ... end block, used to
    specify the configuration of a system of stochastic differential
    equations (SDEs), whose invariant is the Dirichlet distribution. For more
    details on the Dirichlet SDE, see http://dx.doi.org/10.1155/2013/842981.
    Keywords allowed in a dirichlet ... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_b::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappa::string() + "\'. "
    + R"(For an example dirichlet ... end block, see
      doc/html/walker_example_dirichlet.html.)";
  }
};
using dirichlet = keyword< dirichlet_info,  d,i,r,i,c,h,l,e,t >;

struct gendir_info {
  static std::string name() { return "Generalized Dirichlet"; }
  static std::string shortDescription() { return
    "Start configuration block for the generalized Dirichlet SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the gendir ... end
    block, used to specify the configuration of a system of stochastic
    differential equations (SDEs), whose invariant is Lochner's generalized
    Dirichlet distribution. For more details on the generalized Dirichlet
    SDE, see http://dx.doi.org/10.1063/1.4822416. Keywords allowed in a gendir
    ... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_b::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_c::string() + "\', \'"
    + sde_kappa::string() + "\'. "
    + R"(For an example gendir... end block, see
      doc/html/walker_example_gendir.html.)";
  }
};
using gendir = keyword< gendir_info, g,e,n,d,i,r >;

struct wrightfisher_info {
  static std::string name() { return "Wright-Fisher"; }
  static std::string shortDescription() { return
    "Start configuration block for the Wright-Fisher SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the wright_fisher ... end block, used
    to specify the configuration of a system of stochastic differential
    equations (SDEs), whose invariant is the Dirichlet distribution. For more
    details on the Wright-Fisher SDE, see
    http://www.sciencedirect.com/science/article/pii/S0040580912001013.
    Keywords allowed in a wright-fisher... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_omega::string() + "\'. "
    + R"(For an example wright-fisher ... end block, see
      doc/html/walker_example_wf.html.)";
  }
};
using wrightfisher =
  keyword< wrightfisher_info, w,r,i,g,h,t,'-',f,i,s,h,e,r >;

struct skewnormal_info {
  static std::string name() { return "Skew-Normal"; }
  static std::string shortDescription() { return
    "Start configuration block for the Skew-normal SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the skew-normal ... end block, used
    to specify the configuration of a system of stochastic differential
    equations (SDEs), whose invariant is the joint skew-normal distribution.
    For more details on the skew-normal distribution, see
    http://www.jstor.org/stable/2337278. Keywords allowed in an skew-normal ...
    end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_sigmasq::string() + "\', \'"
    + sde_T::string() + "\', \'"
    + sde_lambda::string() + "\'. "
    + R"(For an example skew-normal... end block, see
      doc/html/walker_example_skewnormal.html.)";
  }
};
using skewnormal = keyword< skewnormal_info,  s,k,e,w,'-',n,o,r,m,a,l >;

struct beta_info {
  static std::string name() { return "Beta"; }
  static std::string shortDescription() { return
    "Introduce the beta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the beta ... end block, used to specify
    the configuration of a system of stochastic differential equations
    (SDEs), with linear drift and quadratic diagonal diffusion, whose
    invariant is the joint beta distribution. For more details
    on the beta SDE, see http://doi.org/10.1080/14685248.2010.510843. Keywords
    allowed in a beta ... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_b::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappa::string() + "\'. "
    + R"(For an example beta ... end block, see
      doc/html/walker_example_beta.html.)";
  }
};
using beta = keyword< beta_info, b,e,t,a >;

struct gamma_info {
  static std::string name() { return "Gamma"; }
  static std::string shortDescription() { return
    "Introduce the gamma SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the gamma ... end block, used to
    specify the configuration of a system of stochastic differential equations
    (SDEs), with linear drift and linear diagonal diffusion, whose invariant
    is the joint gamma distribution.)";
  }
};
using gamma = keyword< gamma_info, g,a,m,m,a >;

struct ornstein_uhlenbeck_info {
  static std::string name() { return "Ornstein-Uhlenbeck"; }
  static std::string shortDescription() { return
    "Introduce the Ornstein-Uhlenbeck SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the ornstein-uhlenbeck ... end block,
    used to specify the configuration of a system of stochastic differential
    equations (SDEs), with linear drift and constant diffusion, whose
    invariant is the joint normal distribution. Keywords allowed in an
    ornstein-uhlenbeck ... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_sigmasq::string() + "\', \'"
    + sde_theta::string() + "\', \'"
    + sde_mu::string() + "\'. "
    + R"(For an example ornstein-uhlenbeck ... end block, see
      doc/html/walker_example_ou.html.)";
  }
};
using ornstein_uhlenbeck =
  keyword< ornstein_uhlenbeck_info, o,r,n,s,t,e,i,n,'-',u,h,l,e,n,b,e,c,k >;

struct diag_ou_info {
  static std::string name() { return "Diagonal Ornstein-Uhlenbeck"; }
  static std::string shortDescription() { return
    "Introduce the diagonal Ornstein-Uhlenbeck SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the diag_ou ... end
    block, where 'diag_ou' stands for diagonal Ornstein-Uhlenbeck' and is used
    to specify the configuration of a system of stochastic differential
    equations (SDEs), with linear drift and constant diagonal diffusion, whose
    invariant is the joint normal distribution. Keywords
    allowed in a diagou ... end block: )" + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_sigmasq::string() + "\', \'"
    + sde_theta::string() + "\', \'"
    + sde_mu::string() + "\'. "
    + R"(For an example diagou ... end block, see
      doc/html/walker_example_diagou.html.)";
  }
};
using diag_ou = keyword< diag_ou_info, d,i,a,g,'_',o,u >;

struct control_info {
  static std::string name() { return "control"; }
  static std::string shortDescription()
  { return "Specify the control file name [REQUIRED]"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the control file from which
    detailed user input is parsed.)";
  }
  using alias = Alias< c >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using control = keyword< control_info, c,o,n,t,r,o,l >;

struct smallcrush_info {
  static std::string name() { return "SmallCrush"; }
  static std::string shortDescription() {
    return "Select RNG battery SmallCrush"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the description of the random number
    generator test suite, i.e., battery, 'SmallCrush'. SmallCrush is a
    battery of relatively small number, O(10), of tests, defined in TestU01,
    a library for the empirical testing of random number generators. For more "
    info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.)";
  }
};
using smallcrush = keyword< smallcrush_info, s,m,a,l,l,c,r,u,s,h >;

struct crush_info {
  static std::string name() { return "Crush"; }
  static std::string shortDescription() { return
    "Select RNG battery Crush"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the description of the random number
    generator test suite, i.e., battery, 'Crush'. Crush is a suite of
    stringent statistical tests, O(100), defined in TestU01, a library for
    the empirical testing of random number generators. For more info, see
    http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.)";
  }
};
using crush = keyword< crush_info, c,r,u,s,h >;

struct bigcrush_info {
  static std::string name() { return "BigCrush"; }
  static std::string shortDescription() { return
    "Select RNG battery BigCrush"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the description of the random number
    generator test suite, i.e., battery, 'BigCrush'. BigCrush is a
    suite of very stringent statistical tests, O(100), defined in TestU01, a
    library for the empirical testing of random number generators. For more
    info, see http://www.iro.umontreal.ca/~simardr/testu01/tu01.html.)";
  }
};
using bigcrush = keyword< bigcrush_info, b,i,g,c,r,u,s,h >;

struct verbose_info {
  static std::string name() { return "verbose"; }
  static std::string shortDescription() { return
    "Select verbose screen output"; }
  static std::string longDescription() { return
    R"(This keyword is used to select verbose screen-output as opposed to the
    default quiet output. With quiet output only the most important messages
    are echoed to screen.)";
  }
  using alias = Alias< v >;
};
using verbose = keyword< verbose_info, v,e,r,b,o,s,e >;

struct virtualization_info {
  static std::string name() { return "virtualization"; }
  static std::string shortDescription() { return
    R"(Set degree of virtualization)"; }
  static std::string longDescription() { return
    R"(This option is used to set the degree of virtualization
    (over-decomposition). The virtualization parameter, is a real number
    between 0.0 and 1.0, inclusive, which controls the degree of
    virtualization or over-decomposition. Independent of the value of
    virtualization the work is approximately evenly distributed among the
    available processing elements. For zero virtualization (no
    over-decomposition), the work is simply decomposed into
    total_work/numPEs, which yields the smallest number of Charm++ chares and
    the largest chunks of work units. The other extreme is unity
    virtualization, which decomposes the total work into the smallest size
    work units possible, yielding the largest number of Charm++ chares.
    Obviously, the optimum will be between 0.0 and 1.0, depending on the
    problem.)";
  }
  using alias = Alias< u >;
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static constexpr type upper = 1.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "real between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using virtualization =
  keyword< virtualization_info, v,i,r,t,u,a,l,i,z,a,t,i,o,n >;

struct pdf_info {
  static std::string name() { return "pdf"; }
  static std::string shortDescription() { return
    "Specify the name of the PDF output file"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the output file in which to
    store probability density functions (PDFs) during a simulation.)";
  }
  using alias = Alias< p >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using pdf = keyword< pdf_info, p,d,f >;

struct stat_info {
  static std::string name() { return "stat"; }
  static std::string shortDescription() { return
    "Specify the name of the statistical moments output file"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the name of the output file in which to
    store statistical moments during a simulation.)";
  }
  using alias = Alias< s >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using stat = keyword< stat_info, s,t,a,t >;

struct input_info {
  static std::string name() { return "input"; }
  static std::string shortDescription() { return "Specify the input file"; }
  static std::string longDescription() { return
    R"(This option is used to define the input file name.)";
  }
  using alias = Alias< i >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using input = keyword< input_info, i,n,p,u,t >;

struct output_info {
  static std::string name() { return "output"; }
  static std::string shortDescription() { return "Specify the output file"; }
  static std::string longDescription() { return
    R"(This option is used to define the output file name.)";
  }
  using alias = Alias< o >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using output = keyword< output_info, o,u,t,p,u,t >;

////////// USED ONLY BY EXECUTABLE 'quinoa', NOT YET FULLY DOCUMENTED //////////

struct mix_iem_info {
  static std::string name() { return "IEM"; }
  static std::string shortDescription() { return
    "Interaction by exchange with the mean"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_iem', is short for interaction by exchange with
    the mean (IEM). It is a relaxation-type material mix model intended
    shear-driven flows.)";
  }
};
using mix_iem = keyword<mix_iem_info,  m,i,x,'_',i,e,m >;

struct mix_iecm_info {
  static std::string name() { return "IECM"; }
  static std::string shortDescription()
  { return "Interaction by exchange with the conditional mean"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_iecm', is short for interaction by exchange with
    the conditional mean (IECM). It is a relaxation-type material mix model
    intended shear-driven flows.)";
  }
};
using mix_iecm = keyword<mix_iecm_info,  m,i,x,'_',i,e,c,m >;

struct mix_dir_info {
  static std::string name() { return "Dirichlet"; }
  static std::string shortDescription() { return "Dirichlet"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_dir', is short for Dirichlet. It is a material
    mix model that explicitly satisfies the unit-sum requirement for all
    statistical samples.)";
  }
};
using mix_dir = keyword<mix_dir_info,  m,i,x,'_',d,i,r >;

struct mix_gendir_info {
  static std::string name() { return "Generalized Dirichlet"; }
  static std::string shortDescription() { return "Generalized Dirichlet"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_gendir', is short for Lochner's generalized
    Dirichlet. It is a material mix model that explicitly satisfies the
    unit-sum requirement for all statistical samples.)";
  }
};
using mix_gendir = keyword<mix_gendir_info,  m,i,x,'_',g,e,n,d,i,r >;

struct hommix_info {
  static std::string name() { return "HomMix"; }
  static std::string shortDescription()
  { return "Homogeneous material mixing"; }
  static std::string longDescription() { return
    R"(Physics option, 'hommix', is short for homogeneous material mixing. It is
    the simplest physics option that can be used to research, develop, and
    test material mixing models independent, i.e., decoupled from other
    equations. Only a set of scalar equations are advanced which can be
    coupled to each other. The keyword 'hommix' introduces the hommix ... end
    block, selecting and describing the parameters of the mixing model(s).
    The common physics keywords are recognized.)";
  }
};
using hommix = keyword<hommix_info, h,o,m,m,i,x>;

struct homhydro_info {
  static std::string name() { return "HomHydro"; }
  static std::string shortDescription() { return "Homogeneous hydrodynamics"; }
  static std::string longDescription() { return
    R"(Physics option, 'homhydro', is short for homogeneous hydrodynamics. It is
    the simplest physics option that can be used to research, develop, and
    test hydrodynamics models independent of, i.e., decoupled from other
    equations. Only a set of momentum equations are advanced whose components
    can be coupled to each other. The keyword 'homhydro' introduces the
    homhydro ... end block, selecting and describing the parameters of the
    hydrodynamics model(s). The common physics keywords are recognized.)";
  }
};
using homhydro = keyword<homhydro_info,  h,o,m,h,y,d,r,o >;

struct homrt_info {
  static std::string name() { return "HomRT"; }
  static std::string shortDescription()
  { return "Homogeneous Rayleigh-Taylor"; }
  static std::string longDescription() { return
    R"(Physics option, 'homrt', is short for homogeneous Rayleigh-Taylor. It is
    the simplest physics option that can be used to research, develop, and
    test hydrodynamics models for variable-density hydrodynamics and coupled
    material mixing, independent, i.e., decoupled from other equations. Only
    a set of mass and momentum conservation equations are advanced whose
    components can be coupled to each other. The keyword 'homrt' introduces
    the homrt ... end block, selecting and describing the parameters of the
    mass and hydrodynamics model(s). The common physics keywords are
    recognized.)";
  }
};
using homrt = keyword<homrt_info,  h,o,m,r,t >;

struct spinsflow_info {
  static std::string name() { return "SPINSFlow"; }
  static std::string shortDescription() { return
    "Standalone-particle incompressible Navier-Stokes flow";
  }
  static std::string longDescription() { return
    R"(Physics option, 'spinsflow', is short for standalone-particle
    incompressible Navier-Stokes flow. It is a physics option intended for
    inhomogeneous constant-density flow. The transport equations solved are
    the momentum and optionally, energy, and a set of scalars. The
    divergence-constraint is enforced by solving a Poisson equation and
    projection scheme. The keyword 'spinsflow' introduces the spinsflow ...
    end block, selecting and describing the parameters of the above transport
    equations and their models. The common physics keywords are recognized.)";
  }
};
using spinsflow = keyword<spinsflow_info,  s,p,i,n,s,f,l,o,w >;

// This will go away once all the keywords below are documented
struct undefined_info {
  static std::string name() { return "undef"; }
  static std::string shortDescription() { return "undefined"; }
  static std::string longDescription() { return "Undefined."; }
};

using dmpi = keyword<undefined_info,  d,m,p,i >;
using glbi = keyword<undefined_info,  g,l,b,i >;
using glob = keyword<undefined_info, g,l,o,b >;
using pos_inviscid = keyword<undefined_info,  p,o,s,'_',i,n,v,i,s,c,i,d >;
using pos_viscous = keyword<undefined_info,  p,o,s,'_',v,i,s,c,o,u,s >;
using mass_beta = keyword<undefined_info,  m,a,s,s,'_',b,e,t,a >;
using hydro_slm = keyword<undefined_info,  h,y,d,r,o,'_',s,l,m >;
using hydro_glm = keyword<undefined_info,  h,y,d,r,o,'_',g,l,m >;
using mixrate_gamma = keyword<undefined_info,  m,i,x,r,a,t,e,'_',g,a,m,m,a >;
using freq_gamma = keyword<undefined_info,  f,r,e,q,'_',g,a,m,m,a >;
using position = keyword<undefined_info,  p,o,s,i,t,i,o,n >;
using hydro = keyword<undefined_info,  h,y,d,r,o >;
using mix = keyword<undefined_info,  m,i,x >;
using nposition = keyword<undefined_info,  n,p,o,s,i,t,i,o,n >;
using ndensity = keyword<undefined_info,  n,d,e,n,s,i,t,y >;
using nvelocity = keyword<undefined_info,  n,v,e,l,o,c,i,t,y >;
using nscalar = keyword<undefined_info,  n,s,c,a,l,a,r >;
using nfreq = keyword<undefined_info,  n,f,r,e,q >;
using SLM_C0 = keyword<undefined_info,  C,'0' >;
using freq_gamma_C1 = keyword<undefined_info,  C,'1' >;
using freq_gamma_C2 = keyword<undefined_info,  C,'2' >;
using freq_gamma_C3 = keyword<undefined_info,  C,'3' >;
using freq_gamma_C4 = keyword<undefined_info,  C,'4' >;
using Beta_At = keyword<undefined_info,  A,t >;
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

} // kw::

#endif // Keywords_h
