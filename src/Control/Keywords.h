// *****************************************************************************
/*!
  \file      src/Control/Keywords.h
  \author    J. Bakosi
  \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
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
    and long descriptions, defined by member functions _name()_,
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
        // <tpl_install_dir>/include/pegtl/pegtl/constants.hh for examples
        // of what can be passed as template arguments (basically single
        // characters). The one below defines the character 'c' as the alias.
        // Aliases are single character long. The idea of an alias is to have a
        // long as well as a short keyword for the same functionality.
        // Currently, this is only hooked up for command-line arguments and not
        // for control-file keywords, which is intentional. Command-line
        // arguments are a lot less than control file keywords and are more
        // frequently typed by the user. Thus command-line argument aliases are
        // user-friendly. There are many control file keywords and aliases would
        // only cause confusion. Defining an alias for a command-line argument
        // enables the command-line parser to match on '--longer_keyword' as
        // well as on '-c'. Depending on whether the alias typedef is defined
        // for a keyword or not, the correct grammar is automatically generated
        // at compile-time, matching on both the longer keyword as well as on
        // the alias. Defining an alias for a control file keyword can be done
        // but has no effect in a control file parser.
        using alias = Alias< c >;

        // Optional single-character (policy) code. See also kw::Code and
        // <tpl_install_dir/include/pegtl/pegtl/constants.hh for examples
        // of what can be passed as template arguments (basically single
        // characters). The one below defines the character 'C' as the (policy)
        // code. This code is used for abbreviating policies used to configure
        // various orthogonal behaviors of classes using policy-based design.
        using code = Code< C >;

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
// *****************************************************************************
#ifndef Keywords_h
#define Keywords_h

#include <pegtl/contrib/alphabet.hh>

#include "Types.h"
#include "Keyword.h"
#include "QuinoaConfig.h"

//! Keywords used by all input deck and command line parsers
namespace kw {

using namespace pegtl::alphabet;

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
using title = keyword< title_info, pegtl_string_t("title") >;

struct end_info {
  static std::string name() { return "end"; }
  static std::string shortDescription() { return "End of an input block"; }
  static std::string longDescription() { return
    R"(The end of a block is given by the 'end' keyword in the input file.
    Example: "rngs ... end".)";
  }
};
using end = keyword< end_info, pegtl_string_t("end") >;

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
using help = keyword< help_info, pegtl_string_t("help") >;

struct helpctr_info {
  static std::string name() { return "helpctr"; }
  static std::string shortDescription() { return
    "Display one-liner help on all control file keywords"; }
  static std::string longDescription() { return
    R"(This keyword can be used to get a short one-liner help on all control
    file keywords from an executable.)";
  }
  using alias = Alias< C >;
};
using helpctr = keyword< helpctr_info, pegtl_string_t("helpctr") >;

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
using helpkw = keyword< helpkw_info, pegtl_string_t("helpkw") >;

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
using seed = keyword< seed_info, pegtl_string_t("seed") >;

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
using mkl_mcg31 = keyword< mkl_mcg31_info, pegtl_string_t("mkl_mcg31") >;

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
using mkl_r250 = keyword< mkl_r250_info, pegtl_string_t("mkl_r250") >;

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
  keyword< mkl_mrg32k3a_info, pegtl_string_t("mkl_mrg32k3a") >;

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

using mkl_mcg59 = keyword< mkl_mcg59_info, pegtl_string_t("mkl_mcg59") >;

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
using mkl_wh = keyword< mkl_wh_info, pegtl_string_t("mkl_wh") >;

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
  keyword< mkl_mt19937_info, pegtl_string_t("mkl_mt19937") >;

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
using mkl_mt2203 = keyword< mkl_mt2203_info, pegtl_string_t("mkl_mt2203") >;

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
  keyword< mkl_sfmt19937_info, pegtl_string_t("mkl_sfmt19937") >;

struct mkl_sobol_info {
  static std::string name() { return "MKL SOBOL"; }
  static std::string shortDescription() { return
    "Select Intel MKL SOBOL RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_SOBOL', a 32-bit Gray code-based
    random number generator, producing low-discrepancy sequences for
    dimensions 1 .le. s .le. 40 with available user-defined dimensions, provided
    by Intel's Math Kernel Library (MKL). For more info on MKL see https://
    software.intel.com/en-us/articles/intel-math-kernel-library-documentation.)";
  }
};
using mkl_sobol = keyword< mkl_sobol_info, pegtl_string_t("mkl_sobol") >;

struct mkl_niederr_info {
  static std::string name() { return "MKL NIEDERR"; }
  static std::string shortDescription() { return
   "Select Intel MKL NIEDERR RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select 'VSL_BRNG_NIEDERR', a 32-bit Gray
    code-based random number generator, producing low-discrepancy sequences
    for dimensions 1 .le. s .le. 318 with available user-defined dimensions,
    provided by Intel's Math Kernel Library (MKL). For more info on MKL see
    https://software.intel.com/en-us/articles/intel-math-kernel-library-
    documentation.)";
  }
};
using mkl_niederr = keyword< mkl_niederr_info, pegtl_string_t("mkl_niederr") >;

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
  keyword< mkl_iabstract_info, pegtl_string_t("mkl_iabstract") >;

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
  keyword< mkl_dabstract_info, pegtl_string_t("mkl_dabstract") >;

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
  keyword< mkl_sabstract_info, pegtl_string_t("mkl_sabstract") >;

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
  keyword< mkl_nondeterm_info, pegtl_string_t("mkl_nondeterm") >;

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
using standard = keyword< standard_info, pegtl_string_t("standard") >;

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
using accurate = keyword< accurate_info, pegtl_string_t("accurate") >;

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
  keyword< uniform_method_info, pegtl_string_t("uniform_method") >;

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
using boxmuller = keyword< boxmuller_info, pegtl_string_t("boxmuller") >;

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
using boxmuller2 = keyword< boxmuller2_info, pegtl_string_t("boxmuller2") >;

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
using icdf = keyword< icdf_info, pegtl_string_t("icdf") >;

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
  keyword< gaussian_method_info, pegtl_string_t("gaussian_method") >;

struct cja_info {
  static std::string name() { return "CJA"; }
  static std::string shortDescription() { return
   "Select the Cheng, Johnk, Atkinson algorithm for sampling a beta"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Cheng-Johnk-Atkinson method used to
    generate beta random numbers using the Intel Math Kernel Library (MKL)
    random number generators. For more info on MKL see
    https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using cja = keyword< cja_info, pegtl_string_t("cja") >;

struct cja_accurate_info {
  static std::string name() { return "CJA accurate"; }
  static std::string shortDescription() { return
   "Select the accurate Cheng, Johnk, Atkinson algorithm for sampling a beta"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the accurate version of the
    Cheng-Johnk-Atkinson method used to generate beta random numbers using the
    Intel Math Kernel Library (MKL) random number generators. For more info on
    MKL see https://software.intel.com/en-us/
    articles/intel-math-kernel-library-documentation.)";
  }
};
using cja_accurate = keyword< cja_accurate_info, pegtl_string_t("cja_accurate") >;

struct beta_method_info {
  static std::string name() { return "Beta method"; }
  static std::string shortDescription() { return
    "Select an Intel MKL beta RNG method"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the method used to generate beta
    random numbers using the Intel Math Kernel Library (MKL) random number
    generators. Valid options are 'cja' and 'cja_accurate'. For
    more info on MKL see https://software.intel.com/en-us/articles/intel-math-
    kernel-library-documentation.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + cja::string() + "\' | \'"
                  + cja_accurate::string() + '\'';
    }
  };
};
using beta_method =
  keyword< beta_method_info, pegtl_string_t("beta_method") >;

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
using rngsse_gm19 = keyword< rngsse_gm19_info, pegtl_string_t("rngsse_gm19") >;

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
using rngsse_gm29 = keyword< rngsse_gm29_info, pegtl_string_t("rngsse_gm29") >;

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
using rngsse_gm31 = keyword< rngsse_gm31_info, pegtl_string_t("rngsse_gm31") >;

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
using rngsse_gm55 = keyword< rngsse_gm55_info, pegtl_string_t("rngsse_gm55") >;

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
using rngsse_gm61 = keyword< rngsse_gm61_info, pegtl_string_t("rngsse_gm61") >;

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
  keyword< rngsse_gq581_info, pegtl_string_t("rngsse_gq58.1") >;

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
  keyword< rngsse_gq583_info, pegtl_string_t("rngsse_gq58.3") >;

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
  keyword< rngsse_gq584_info, pegtl_string_t("rngsse_gq58.4") >;

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
  keyword< rngsse_mt19937_info, pegtl_string_t("rngsse_mt19937") >;

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
  keyword< rngsse_lfsr113_info, pegtl_string_t("rngsse_lfsr113") >;

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
  keyword< rngsse_mrg32k3a_info, pegtl_string_t("rngsse_mrg32k3a") >;

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
using seq_short = keyword< seq_short_info, pegtl_string_t("short") >;

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
using seq_med = keyword< seq_med_info, pegtl_string_t("medium") >;

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
using seq_long = keyword< seq_long_info, pegtl_string_t("long") >;

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
using seqlen = keyword< seqlen_info, pegtl_string_t("seqlen") >;

struct r123_threefry_info {
  static std::string name() { return "Random123 ThreeFry"; }
  static std::string shortDescription() { return
    "Select Random123 ThreeFry RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the ThreeFry generator, related to the
    Threefish block cipher from Skein Hash Function, provided by the Random123
    random number generator library. For more info on Random123 see
    http://dl.acm.org/citation.cfm?doid=2063405.)";
  }
};
using r123_threefry =
  keyword< r123_threefry_info, pegtl_string_t("r123_threefry") >;

struct r123_philox_info {
  static std::string name() { return "Random123 Philox"; }
  static std::string shortDescription() { return
    "Select Random123 Philox RNG"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Philox generator, based on Feistel
    network and integer multiplication, provided by the Random123
    random number generator library. For more info on Random123 see
    http://dl.acm.org/citation.cfm?doid=2063405.)";
  }
};
using r123_philox =
  keyword< r123_philox_info, pegtl_string_t("r123_philox") >;

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
using pdfs = keyword< pdfs_info, pegtl_string_t("pdfs") >;

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
using txt = keyword< txt_info, pegtl_string_t("txt") >;

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
using gmshtxt = keyword< gmshtxt_info, pegtl_string_t("gmshtxt") >;

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
using gmshbin = keyword< gmshbin_info, pegtl_string_t("gmshbin") >;

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
using exodusii = keyword< exodusii_info, pegtl_string_t("exodusii") >;

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
using pdf_filetype = keyword< filetype_info, pegtl_string_t("filetype") >;

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
using overwrite = keyword< overwrite_info, pegtl_string_t("overwrite") >;

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
using multiple = keyword< multiple_info, pegtl_string_t("multiple") >;

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
using evolution = keyword< evolution_info, pegtl_string_t("evolution") >;

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
using pdf_policy = keyword< policy_info, pegtl_string_t("policy") >;

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
using txt_float_default = keyword< txt_float_default_info, pegtl_string_t("default") >;

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
  keyword< txt_float_scientific_info, pegtl_string_t("scientific") >;

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
using txt_float_fixed = keyword< txt_float_fixed_info, pegtl_string_t("fixed") >;

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
using txt_float_format = keyword< txt_float_format_info, pegtl_string_t("format") >;

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
using precision = keyword< precision_info, pegtl_string_t("precision") >;

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
using elem = keyword< elem_info, pegtl_string_t("elem") >;

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
using node = keyword< node_info, pegtl_string_t("node") >;

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
using pdf_centering = keyword< centering_info, pegtl_string_t("centering") >;

struct raw_info {
  using code = Code< R >;
  static std::string name() { return "raw"; }
  static std::string shortDescription() { return
    "Select the raw initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the raw
    initialization policy. The initialization policy is used to specify how
    the initial conditions are set at t = 0 before time-integration.
    Example: "init raw", which selects raw initialization policy, which
    leaves the memory uninitialized. Note that this option may behave
    differently depending on the particular equation or physical model. See the
    the init policies in DiffEq/InitPolicy.h for valid options.)"; }
};
using raw = keyword< raw_info, pegtl_string_t("raw") >;

struct zero_info {
  using code = Code< Z >;
  static std::string name() { return "zero"; }
  static std::string shortDescription() { return
    "Select the zero initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the zero
    initialization policy. The initialization policy is used to specify how
    the initial conditions are set at t = 0 before time-integration.
    Example: "init zero", which selects zero initialization policy, which
    puts zeros in memory. Note that this option may behave differently
    depending on the particular equation or physical model. See the init
    policies in DiffEq/InitPolicy.h for valid options.)"; }
};
using zero = keyword< zero_info, pegtl_string_t("zero") >;

struct jointdelta_info {
  using code = Code< D >;
  static std::string name() { return "delta"; }
  static std::string shortDescription() { return
    "Select the joint delta initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the joint delta initialization policy.
    The initialization policy is used to specify how the initial conditions are
    set at t = 0 before time-integration. Example: "init zero", which selects
    zero initialization policy, which puts zeros in memory. Note that this
    option may behave differently depending on the particular equation or
    physical model. For an example, see tk::InitPolicies in
    DiffEq/InitPolicy.h for valid options.) The joint delta initialization
    policy can be used to prescribe delta-spikes on the sample space with given
    heights, i.e., probabilities. Example: "init jointdelta" - select delta
    init-policy, "delta spike 0.1 0.3 0.8 0.7 end end" - prescribe two
    delta-spikes at sample space positions 0.1 and 0.8 with spike heights 0.3
    and 0.7, respectively. Note that the sum of the heights must add up to
    unity. See also the help on keyword spike.)"; }
};
using jointdelta = keyword< jointdelta_info, pegtl_string_t("jointdelta") >;

struct jointbeta_info {
  using code = Code< B >;
  static std::string name() { return "beta"; }
  static std::string shortDescription() { return
    "Select the joint beta initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the joint beta initialization policy.
    The initialization policy is used to specify how the initial conditions are
    set at t = 0 before time-integration. Example: "init zero", which selects
    zero initialization policy, which puts zeros in memory. Note that this
    option may behave differently depending on the particular equation or
    physical model. For an example, see tk::InitPolicies in
    DiffEq/InitPolicy.h for valid options.) The joint beta initialization
    policy can be used to prescribe a multi-dimensional sample space where the
    samples are generated from a joint beta distribution with independent
    marginal univariate beta distributions.)";
  }
};
using jointbeta = keyword< jointbeta_info, pegtl_string_t("jointbeta") >;

struct init_info {
  using code = Code< i >;
  static std::string name() { return "initialization policy"; }
  static std::string shortDescription() { return
    "Select initialization policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select an
    initialization policy. This is used to specify how the initial conditions
    are set at t = 0 before time-integration. Example: "init raw", which
    selects raw initialization policy, which leaves the memory uninitialized.
    Note that this option may behave differently depending on the particular
    equation or physical model. See the init policies in
    DiffEq/InitPolicy.h for valid options.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + raw::string() + "\' | \'"
                  + zero::string() + "\' | \'"
                  + jointdelta::string() + "\' | \'"
                  + jointbeta::string() + '\'';
    }
  };
};
using init = keyword< init_info, pegtl_string_t("init") >;

struct const_info {
  using code = Code< C >;
  static std::string name() { return "constant"; }
  static std::string shortDescription() { return
    "Select constant coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the
    constant coefficients policy. The coefficients policy is used to specify
    how the coefficients are set at each time step during time-integration.
    Example: "coeff const", which selects constant coefficients policy,
    which sets constant coefficients before t = 0 and leaves the coefficients
    unchanged during time integration. Note that this option may behave
    differently depending on the particular equation or physical model.)"; }
};
using constant = keyword< const_info, pegtl_string_t("const") >;

struct decay_info {
  using code = Code< D >;
  static std::string name() { return "decay"; }
  static std::string shortDescription() { return
    "Select decay coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the decay coefficients policy. This policy
    (or model) is used to constrain a beta stochastic differential equation so
    that its variance, <y^2>, always decays. A coefficients policy, in general,
    is used to specify how the coefficients are set at each time step during
    time-integration. Example: "coeff const", which selects constant
    coefficients policy, which sets constant coefficients before t = 0 and
    leaves the coefficients unchanged during time integration. Note that this
    option may behave differently depending on the particular equation or
    physical model.)"; }
};
using decay = keyword< decay_info, pegtl_string_t("decay") >;

struct homdecay_info {
  using code = Code< H >;
  static std::string name() { return "homogeneous decay"; }
  static std::string shortDescription() { return
    "Select homogeneous decay coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the homogeneous decay coefficients policy.
    This policy (or model) is used to constrain a beta stochastic differential
    equation so that its variance, <y^2>, always decays and its mean, <R> =
    rho2/(1+r<RY>/<R>), where Y = <Y> + y, does not change in time. Note that
    R = rho2/(1+rY). This policy is similar to 'montecarlo_homdecay', but
    computes the SDE coefficient S in a different but statistically equivalent
    way. While 'homdecay' only requires the estimation of statistics, <R>,
    <r^2>, and <r^3>, 'montecarlo_homdecay' requires <R^2>, <YR^2>, and
    <Y(1-Y)R^3>. A coefficients policy, in general, is used to specify how the
    coefficients are set at each time step during time-integration. Example:
    "coeff const", which selects constant coefficients policy, which sets
    constant coefficients before t = 0 and leaves the coefficients unchanged
    during time integration. Note that this option may behave differently
    depending on the particular equation or physical model.)"; }
};
using homdecay = keyword< homdecay_info, pegtl_string_t("homdecay") >;

struct montecarlo_homdecay_info {
  using code = Code< M >;
  static std::string name() { return "Monte Carlo homogeneous decay"; }
  static std::string shortDescription() { return
    "Select Monte Carlo homogeneous decay coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Monte Carlo homogeneous decay
    coefficients policy. This policy (or model) is used to constrain a beta
    stochastic differential equation (SDE) so that its variance, <y^2>, always
    decays and its mean, <R> = rho2/(1+r<RY>/<R>), where Y = <Y> + y, does not
    change in time. Note that R = rho2/(1+rY). This policy is similar to
    'homdecay', but computes the the SDE coefficient S in a different but
    statistically equivalent way. While 'homdecay' only requires the estimation
    of statistics, <R>, <r^2>, and <r^3>, 'montecarlo_homdecay' requires <R^2>,
    <YR^2>, and <Y(1-Y)R^3>. A coefficients policy, in general, is used to
    specify how the coefficients are set at each time step during
    time-integration. Example: "coeff const", which selects constant
    coefficients policy, which sets constant coefficients before t = 0 and
    leaves the coefficients unchanged during time integration. Note that this
    option may behave differently depending on the particular equation or
    physical model.)"; }
};
using montecarlo_homdecay =
  keyword< montecarlo_homdecay_info, pegtl_string_t("montecarlo_homdecay") >;

struct hydrotimescale_info {
  using code = Code< T >;
  static std::string name() { return "hydro-timescale homogeneous decay"; }
  static std::string shortDescription() { return
    "Select hydro-timescale homogeneous decay coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the hydrodynamics-timescale homogeneous
    decay coefficients policy. This policy (or model) is used to constrain a
    beta stochastic differential equation (SDE) so that its variance, <y^2>,
    always decays and its mean, <R> = rho2/(1+r<RY>/<R>), where Y = <Y> + y,
    does not change in time. Note that R = rho2/(1+rY). This policy is similar
    to 'homdecay' as well as 'montecarlo_homdecay', but instead of simply
    constraining b' and kappa' to ensure decay in the evolution of <y^2>, b' and
    kappa' are specified as functions of an externally-specified hydrodynamics
    time scale, as a function of time. This policy is more similar to 'homdecay'
    than to 'montecarlo_homdecay' in that only requires the estimation
    of statistics, <R>, <r^2>, and <r^3>. A coefficients policy, in general, is
    used to specify how the coefficients are set at each time step during
    time-integration. Example: "coeff const", which selects constant
    coefficients policy, which sets constant coefficients before t = 0 and
    leaves the coefficients unchanged during time integration. Note that this
    option may behave differently depending on the particular equation or
    physical model.)"; }
};
using hydrotimescale =
  keyword< hydrotimescale_info, pegtl_string_t("hydrotimescale") >;

struct coeff_info {
  using code = Code< c >;
  static std::string name() { return "coeficients fpolicy"; }
  static std::string shortDescription() { return
    "Select the coefficients policy"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a
    coefficients policy. This is used to specify how the coefficients are set
    at each time step during time-integration. Example: "coeff const",
    which selects constant coefficients policy, which sets constant
    coefficients before t = 0 and leaves the coefficients unchanged during
    time integration. Note that this option may behave differently depending
    on the particular equation or physical model.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + constant::string() + '\'';
    }
  };
};
using coeff = keyword< coeff_info,  pegtl_string_t("coeff") >;

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
using walker = keyword< walker_info, pegtl_string_t("walker") >;

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
using npar = keyword< npar_info, pegtl_string_t("npar") >;

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
using nstep = keyword< nstep_info, pegtl_string_t("nstep") >;

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
using term = keyword< term_info, pegtl_string_t("term") >;

struct t0_info {
  static std::string name() { return "t0"; }
  static std::string shortDescription() { return
    "Set starting non-dimensional time"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the starting time in a simulation.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using t0 = keyword< t0_info, pegtl_string_t("t0") >;

struct dt_info {
  static std::string name() { return "dt"; }
  static std::string shortDescription() { return
    "Select constant time step size"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the time step size that used as a
    constant during simulation. Setting 'cfl' and 'dt' are mutually
    exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using dt = keyword< dt_info, pegtl_string_t("dt") >;

struct cfl_info {
  static std::string name() { return "CFL"; }
  static std::string shortDescription() { return
    "Set the Courant-Friedrichs-Lewy (CFL) coefficient"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the CFL coefficient for
    variable-time-step-size simulations. Setting 'cfl' and 'dt' are mutually
    exclusive. If both 'cfl' and 'dt' are set, 'dt' wins.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using cfl = keyword< cfl_info, pegtl_string_t("cfl") >;

struct ncomp_info {
  static std::string name() { return "ncomp"; }
  static std::string shortDescription() { return
    "Set number of scalar components for a system of SDEs"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of scalar components of a
    vector.)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 0;
    static std::string description() { return "uint"; }
  };
};
using ncomp = keyword< ncomp_info,  pegtl_string_t("ncomp") >;

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
using ttyi = keyword< ttyi_info, pegtl_string_t("ttyi") >;

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
using interval = keyword< interval_info, pegtl_string_t("interval") >;

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
using statistics = keyword< statistics_info, pegtl_string_t("statistics") >;

struct plotvar_info {
  static std::string name() { return "plotvar"; }
  static std::string shortDescription() { return
    "Start of plotvar input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to start a block in the input file containing the
    list and settings of requested field output.)";
  }
};
using plotvar = keyword< plotvar_info, pegtl_string_t("plotvar") >;

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
using rngs = keyword< rngs_info, pegtl_string_t("rngs") >;

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
      return '\'' + r123_threefry::string()   + "\' | \'"
                  + r123_philox::string()
                  #ifdef HAS_RNGSSE2
                                              + "\' | \'"
                  + rngsse_gm19::string()     + "\' | \'"
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
                  #endif
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
using rng = keyword< rng_info, pegtl_string_t("rng") >;

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
using sde_omega = keyword< sde_omega_info, pegtl_string_t("omega") >;

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
using sde_b = keyword< sde_b_info,  pegtl_string_t("b") >;

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
using sde_S = keyword< sde_S_info,  pegtl_string_t("S") >;

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
using sde_kappa = keyword< sde_kappa_info,  pegtl_string_t("kappa") >;

struct sde_bprime_info {
  static std::string name() { return "bprime"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) bprime)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "bprime 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = sde_b_info::expect::type;
    static std::string description() { return "real(s)"; }
  };
};
using sde_bprime = keyword< sde_bprime_info,  pegtl_string_t("bprime") >;

struct sde_kappaprime_info {
  static std::string name() { return "kappaprime"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) kappaprime)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "kappaprime 5.0 2.0 3.0 end". The length of the vector depends on the
    particular type of SDE system and is controlled by the preceding keyword
    'ncomp'.)"; }
  struct expect {
    using type = sde_kappa_info::expect::type;
    static std::string description() { return "real(s)"; }
  };
};
using sde_kappaprime = keyword< sde_kappaprime_info,  pegtl_string_t("kappaprime") >;

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
using sde_c = keyword< sde_c_info,  pegtl_string_t("c") >;

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
using sde_sigmasq = keyword< sde_sigmasq_info, pegtl_string_t("sigmasq") >;

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
using sde_theta = keyword< sde_theta_info, pegtl_string_t("theta") >;

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
using sde_mu = keyword< sde_mu_info, pegtl_string_t("mu") >;

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
using sde_T = keyword< sde_T_info, pegtl_string_t("T") >;

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
using sde_lambda = keyword< sde_lambda_info, pegtl_string_t("lambda") >;

struct spike_info {
  static std::string name() { return "spike"; }
  static std::string shortDescription() { return
    R"(Configure a delta spike)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the configuration of delta spikes for,
    the delta initialization policy. The configuration is given by an even set
    of real numbers inside a spike...end block. Example: "spike 0.1 1.0 end",
    which specifies a delta spike at sample space position 0.1 with relative
    height 1.0. The height must be between [0.0...1.0] inclusive and specifies a
    relative probability. See also the help on keyword icdelta.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "even reals"; }
  };
};
using spike = keyword< spike_info, pegtl_string_t("spike") >;

struct icdelta_info {
  static std::string name() { return "icdelta"; }
  static std::string shortDescription() { return
    R"(Introduce a icdelta...end block used to configure delta spikes)"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a icdelta...end block in which delta
    spikes are configured for the delta initialization policy. Example:
    "init jointdelta" - select joint delta init-policy,"icdelta spike 0.1 0.3
    0.9 0.7 end end" - prescribe a univariate distribution that consists of two
    delta-spikes at sample space positions 0.1 and 0.9 with spike heights 0.3
    and 0.7, respectively. Note that the sum of the heights must add up to
    unity. See also the help on keyword jointdelta and spike.)"; }
};
using icdelta = keyword< icdelta_info, pegtl_string_t("icdelta") >;

struct betapdf_info {
  static std::string name() { return "betapdf"; }
  static std::string shortDescription() { return
    R"(Configure a beta distribution)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the configuration of beta distributions
    for the beta initialization policy. The configuration is given by four
    real numbers inside a betapdf...end block. Example: "betapdf 0.2 0.3 0.0 1.0
    end", which specifies a univariate beta distribution with shape parameters
    0.2 and 0.3, displacement 0.0, and scale 1.0. See also the help on keyword
    icbeta.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "4 reals"; }
  };
};
using betapdf = keyword< betapdf_info, pegtl_string_t("betapdf") >;

struct icbeta_info {
  static std::string name() { return "icbeta"; }
  static std::string shortDescription() { return
    R"(Introduce an icbeta...end block used to configure beta distributions)"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an icbeta...end block in which beta
    distributions are configured for the beta initialization policy. Example:
    "init jointbeta" - select beta init-policy,"icbeta betapdf 0.2 0.3 0.0 1.0
    end end" - prescribe a univariate beta distribution with shape parameters
    0.2 and 0.3, displacement 0.0, and scale 1.0. See also the help on keyword
    jointbeta and betapdf.)"; }
};
using icbeta = keyword< icbeta_info, pegtl_string_t("icbeta") >;

struct ic_info {
  static std::string name() { return "ic"; }
  static std::string shortDescription() { return
    R"(Introduce an ic...end block used to configure initial conditions)"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an ic...end block used to set initial
    conditions. Example: "ic density 1.0 end" - set the initial density field to
    1.0 across the whole domain.)"; }
};
using ic = keyword< ic_info, pegtl_string_t("i,c") >;

struct depvar_info {
  static std::string name() { return "depvar"; }
  static std::string shortDescription() { return
    "Select dependent variable (in a relevant block)"; }
  static std::string longDescription() { return
    R"(Dependent variable, e.g, in differential equations.)"; }
};
using depvar = keyword< depvar_info, pegtl_string_t("depvar") >;

struct sde_rho2_info {
  static std::string name() { return "rho2"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) rho2)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "rho2 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_rho2 = keyword< sde_rho2_info,  pegtl_string_t("rho2") >;

struct sde_rcomma_info {
  static std::string name() { return "rcomma"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) rcomma)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "rcomma 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_rcomma = keyword< sde_rcomma_info, pegtl_string_t("rcomma") >;

struct sde_r_info {
  static std::string name() { return "r"; }
  static std::string shortDescription() { return
    R"(Set SDE parameter(s) r)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of stochastic differential equations. Example:
    "r 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of SDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using sde_r = keyword< sde_r_info, pegtl_string_t("r") >;

struct dirichlet_info {
  static std::string name() { return "Dirichlet"; }
  static std::string shortDescription() { return
    "Start configuration block for the Dirichlet SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a dirichlet ... end block, used to
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
using dirichlet = keyword< dirichlet_info,  pegtl_string_t("dirichlet") >;

struct gendir_info {
  static std::string name() { return "Generalized Dirichlet"; }
  static std::string shortDescription() { return
    "Start configuration block for the generalized Dirichlet SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a gendir ... end
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
using gendir = keyword< gendir_info, pegtl_string_t("gendir") >;

struct wrightfisher_info {
  static std::string name() { return "Wright-Fisher"; }
  static std::string shortDescription() { return
    "Start configuration block for the Wright-Fisher SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a wright_fisher ... end block, used
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
  keyword< wrightfisher_info, pegtl_string_t("wright-fisher") >;

struct skewnormal_info {
  static std::string name() { return "Skew-Normal"; }
  static std::string shortDescription() { return
    "Start configuration block for the Skew-normal SDE"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a skew-normal ... end block, used
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
using skewnormal = keyword< skewnormal_info,  pegtl_string_t("skew-normal") >;

struct beta_info {
  static std::string name() { return "Beta"; }
  static std::string shortDescription() { return
    "Introduce the beta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a beta ... end block, used to specify
    the configuration of a system of stochastic differential equations (SDEs),
    with linear drift and quadratic diagonal diffusion, whose invariant is the
    joint beta distribution. For more details on the beta SDE, see
    http://doi.org/10.1080/14685248.2010.510843 and src/DiffEq/Beta.h. Keywords
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
using beta = keyword< beta_info, pegtl_string_t("beta") >;

struct numfracbeta_info {
  static std::string name() { return "Number-fraction beta"; }
  static std::string shortDescription() { return
    "Introduce the numfracbeta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a numfracbeta ... end block, used to
    specify the configuration of a system of number-fraction beta SDEs, a system
    of stochastic differential equations (SDEs), in which, in addition to the
    dependent variable, computed with linear drift and quadratic diagonal
    diffusion (whose invariant is joint beta), two additional variables are
    computed. In other words, this is a beta SDE but there are two additional
    stochastic variables computed based on the beta SDE. If X is governed by the
    beta SDE, then the number-fraction beta SDE additionally governs rho(X) and
    V(X), where both rho and V are random variables, computed by rho(X) = rho2
    ( 1 - r' X ), and V(X) = 1 / [ rho2 ( 1 - r'X ) ]. For more details on the
    beta SDE, see http://doi.org/10.1080/14685248.2010.510843 and
    src/DiffEq/Beta.h. Keywords allowed in a numfracbeta ... end block: )"
    + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_b::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappa::string() + "\', \'"
    + sde_rho2::string() + "\', \'"
    + sde_rcomma::string() + "\'. "
    + R"(For an example numfracbeta ... end block, see
      doc/html/walker_example_numfracbeta.html.)";
  }
};
using numfracbeta = keyword< numfracbeta_info, pegtl_string_t("numfracbeta") >;

struct massfracbeta_info {
  static std::string name() { return "Mass-fraction beta"; }
  static std::string shortDescription() { return
    "Introduce the massfracbeta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a massfracbeta ... end block, used to
    specify the configuration of a system of number-fraction beta SDEs, a system
    of stochastic differential equations (SDEs), in which, in addition to the
    dependent variable, computed with linear drift and quadratic diagonal
    diffusion (whose invariant is joint beta), two additional variables are
    computed. In other words, this is a beta SDE but there are two additional
    stochastic variables computed based on the beta SDE. If Y is governed by the
    beta SDE, then the mass-fraction beta SDE additionally governs rho(Y) and
    V(Y), where both rho and V are random variables, computed by rho(Y) = rho2 /
    ( 1 + r Y ), and V(Y) = ( 1 + r Y ) / rho2. For more details on the beta
    SDE, see http://doi.org/10.1080/14685248.2010.510843 and src/DiffEq/Beta.h.
    Keywords allowed in a massfracbeta ... end block: )"
    + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_b::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappa::string() + "\', \'"
    + sde_rho2::string() + "\', \'"
    + sde_r::string() + "\'. "
    + R"(For an example massfracbeta ... end block, see
      doc/html/walker_example_massfracbeta.html.)";
  }
};
using massfracbeta = keyword< massfracbeta_info, pegtl_string_t("massfracbeta") >;

struct mixnumfracbeta_info {
  static std::string name() { return "Mix number-fraction beta"; }
  static std::string shortDescription() { return
    "Introduce the mixnumfracbeta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a mixnumfracbeta ... end block, used
    to specify the configuration of a system of mix number-fraction beta SDEs, a
    system of stochastic differential equations (SDEs), whose solution is the
    joint beta distribution and in which the usual beta SDE parameters b and
    kappa are specified via functions that constrain the beta SDE to be
    consistent with the turbulent mixing process. The mix number-fraction beta
    SDE is similar to the number-fraction beta SDE, only the process is made
    consistent with the no-mix and fully mixed limits via the specification of
    the SDE coefficients b and kappa. As in the number-fraction beta SDE, X is
    governed by the beta SDE and two additional stochastic variables are
    computed. However, in the mix number-fraction beta SDE the parameters b and
    kappa are given by b = Theta * b' and kappa = kappa' * <x^2>, where Theta =
    1 - <x^2> / [ <X> ( 1 - <X> ], the fluctuation about the mean, <X>, is
    defined as usual: x = X - <X>, and b' and kappa' are user-specified
    constants. Similar to the number-fraction beta SDE, there two additional
    random variables computed besides, X, and they are rho(X) and V(X). For more
    detail on the number-fraction beta SDE, see the help on keyword
    'numfracbeta'. For more details on the beta SDE, see
    http://doi.org/10.1080/14685248.2010.510843 and src/DiffEq/Beta.h. Keywords
    allowed in a mixnumfracbeta ... end block: )"
    + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_bprime::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappaprime::string() + "\', \'"
    + sde_rho2::string() + "\', \'"
    + sde_rcomma::string() + "\'. "
    + R"(For an example mixnumfracbeta ... end block, see
      doc/html/walker_example_mixnumfracbeta.html.)";
  }
};
using mixnumfracbeta =
  keyword< mixnumfracbeta_info, pegtl_string_t("mixnumfracbeta") >;

struct mixmassfracbeta_info {
  static std::string name() { return "Mix mass-fraction beta"; }
  static std::string shortDescription() { return
    "Introduce the mixmassfracbeta SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a mixmassfracbeta ... end block, used
    to specify the configuration of a system of mix mass-fraction beta SDEs, a
    system of stochastic differential equations (SDEs), whose solution is the
    joint beta distribution and in which the usual beta SDE parameters b and
    kappa are specified via functions that constrain the beta SDE to be
    consistent with the turbulent mixing process. The mix mass-fraction beta
    SDE is similar to the mass-fraction beta SDE, only the process is made
    consistent with the no-mix and fully mixed limits via the specification of
    the SDE coefficients b and kappa. As in the mass-fraction beta SDE, Y is
    governed by the beta SDE and two additional stochastic variables are
    computed. However, in the mix mass-fraction beta SDE the parameters b and
    kappa are given by b = Theta * b' and kappa = kappa' * <y^2>, where Theta =
    1 - <y^2> / [ <Y> ( 1 - <Y> ], the fluctuation about the mean, <Y>, is
    defined as usual: y = Y - <Y>, and b' and kappa' are user-specified
    constants. Similar to the mass-fraction beta SDE, there two additional
    random variables computed besides, Y, and they are rho(Y) and V(Y). For more
    detail on the mass-fraction beta SDE, see the help on keyword
    'massfracbeta'. For more details on the beta SDE, see
    http://doi.org/10.1080/14685248.2010.510843 and src/DiffEq/Beta.h. Keywords
    allowed in a mixmassfracbeta ... end block: )"
    + std::string("\'")
    + depvar::string()+ "\', \'"
    + ncomp::string() + "\', \'"
    + rng::string() + "\', \'"
    + init::string() + "\', \'"
    + coeff::string() + "\', \'"
    + sde_bprime::string() + "\', \'"
    + sde_S::string() + "\', \'"
    + sde_kappaprime::string() + "\', \'"
    + sde_rho2::string() + "\', \'"
    + sde_r::string() + "\'. "
    + R"(For an example mixmassfracbeta ... end block, see
      doc/html/walker_example_mixmassfracbeta.html.)";
  }
};
using mixmassfracbeta =
  keyword< mixmassfracbeta_info, pegtl_string_t("mixmassfracbeta") >;

struct eq_A005H_info {
  static std::string name() { return "eq_A005H"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.05, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.05, IC: light << heavy.)"; }
};
using eq_A005H = keyword< eq_A005H_info, pegtl_string_t("eq_A005H") >;

struct eq_A005S_info {
  static std::string name() { return "eq_A005S"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.05, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.05, IC: light = heavy.)"; }
};
using eq_A005S = keyword< eq_A005S_info, pegtl_string_t("eq_A005S") >;

struct eq_A005L_info {
  static std::string name() { return "eq_A005L"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.05, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.05, IC: light >> heavy.)"; }
};
using eq_A005L = keyword< eq_A005L_info, pegtl_string_t("eq_A005L") >;

struct eq_A05H_info {
  static std::string name() { return "eq_A05H"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.5, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.5, IC: light << heavy.)"; }
};
using eq_A05H = keyword< eq_A05H_info, pegtl_string_t("eq_A05H") >;

struct eq_A05S_info {
  static std::string name() { return "eq_A05S"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.5, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.5, IC: light = heavy.)"; }
};
using eq_A05S = keyword< eq_A05S_info, pegtl_string_t("eq_A05S") >;

struct eq_A05L_info {
  static std::string name() { return "eq_A05L"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.5, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.5, IC: light >> heavy.)"; }
};
using eq_A05L = keyword< eq_A05L_info, pegtl_string_t("eq_A05L") >;

struct eq_A075H_info {
  static std::string name() { return "eq_A075H"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.75, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.75, IC: light << heavy.)"; }
};
using eq_A075H = keyword< eq_A075H_info, pegtl_string_t("eq_A075H") >;

struct eq_A075S_info {
  static std::string name() { return "eq_A075S"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.75, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.75, IC: light = heavy.)"; }
};
using eq_A075S = keyword< eq_A075S_info, pegtl_string_t("eq_A075S") >;

struct eq_A075L_info {
  static std::string name() { return "eq_A075L"; }
  static std::string shortDescription() { return "Select inverse equilibrium "
   "hydro time scale from DNS of HRT, A=0.75, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Inverse equilibrium hydrodynamics time scale from DNS of homogeneous
       Rayleigh-Taylor instability, tau_eq, A = 0.75, IC: light >> heavy.)"; }
};
using eq_A075L = keyword< eq_A075L_info, pegtl_string_t("eq_A075L") >;

struct hydrotimescales_info {
  static std::string name() { return "hydrotimescales"; }
  static std::string shortDescription() { return
    R"(Set MixMassFractionBeta SDE parameter(s) hydrotimescales)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of strings used to parameterize
    the system of stochastic differential equations, configured in a
    mixmassfracbeta ... end block. Within the mixmassfracbeta ... end block the
    coefficients policy must be set to 'hydrotimescale' in order for the
    hydrotimescales ... end  block to be in effect. The 'hydrotimescales'
    keyword is then used to specify a list of strings, each specifying which
    inverse time scale should be used for the particular component integrated.
    Available time scales are defined in src/DiffEq/HydroTimescales.h. Example:
    "hydrotimescales eq_A05S eq_A05H eq_A05L eq_A05S eq_A05S end", which
    configures five inverse hydrodynamics time scales associated to 5
    components, i.e., 5 scalar stochastic differential equations, integrated,
    specified and configured within the given mixmassfracbeta ... end block. The
    length of the hydrotimescales vector depends on the number of scalar
    components and is controlled by the preceding keyword 'ncomp'. For
    mixmassfracbeta, ncomp is the actual number of scalar components * 4, since
    mixmassfractionbeta always computes 4 additional derived stochastic
    variables (in a diagnostic) fashion. See also MixMassFractionBeta::derived()
    in src/DiffEq/MixMassFractionBeta.h. Keywords allowed in a hydrotimescales
    ... end block: )" + std::string("\'")
    + eq_A005H::string() + "\', \'"
    + eq_A005S::string() + "\', \'"
    + eq_A005L::string() + "\', \'"
    + eq_A05H::string() + "\', \'"
    + eq_A05S::string() + "\', \'"
    + eq_A05L::string() + "\', \'"
    + eq_A075H::string() + "\', \'"
    + eq_A075S::string() + "\', \'"
    + eq_A075L::string() + "\'. "
    + R"(For an example hydrotimescales ... end block, see
      doc/html/walker_example_mixmassfracbeta.html.)"; }
  struct expect {
    using type = std::string;
    static std::string description() { return "string(s)"; }
  };
};
using hydrotimescales =
  keyword< hydrotimescales_info, pegtl_string_t("hydrotimescales") >;

struct prod_A005H_info {
  static std::string name() { return "prod_A005H"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.05, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.05, IC: light << heavy.)"; }
};
using prod_A005H = keyword< prod_A005H_info, pegtl_string_t("prod_A005H") >;

struct prod_A005S_info {
  static std::string name() { return "prod_A005S"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.05, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.05, IC: light = heavy.)"; }
};
using prod_A005S = keyword< prod_A005S_info, pegtl_string_t("prod_A005S") >;

struct prod_A005L_info {
  static std::string name() { return "prod_A005L"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.05, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.05, IC: light >> heavy.)"; }
};
using prod_A005L = keyword< prod_A005L_info, pegtl_string_t("prod_A005L") >;

struct prod_A05H_info {
  static std::string name() { return "prod_A05H"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.5, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.5, IC: light << heavy.)"; }
};
using prod_A05H = keyword< prod_A05H_info, pegtl_string_t("prod_A05H") >;

struct prod_A05S_info {
  static std::string name() { return "prod_A05S"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.5, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.5, IC: light = heavy.)"; }
};
using prod_A05S = keyword< prod_A05S_info, pegtl_string_t("prod_A05S") >;

struct prod_A05L_info {
  static std::string name() { return "prod_A05L"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.5, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.5, IC: light >> heavy.)"; }
};
using prod_A05L = keyword< prod_A05L_info, pegtl_string_t("prod_A05L") >;

struct prod_A075H_info {
  static std::string name() { return "prod_A075H"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.75, IC:light<<heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.75, IC: light << heavy.)"; }
};
using prod_A075H = keyword< prod_A075H_info, pegtl_string_t("prod_A075H") >;

struct prod_A075S_info {
  static std::string name() { return "prod_A075S"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.75, IC:light=heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.75, IC: light = heavy.)"; }
};
using prod_A075S = keyword< prod_A075S_info, pegtl_string_t("prod_A075S") >;

struct prod_A075L_info {
  static std::string name() { return "prod_A075L"; }
  static std::string shortDescription() { return "Select production divided by "
    "dissipation rate from DNS of HRT, A=0.75, IC:light>>heavy"; }
  static std::string longDescription() { return
    R"(Production divided by dissipation rate from DNS of homogeneous
       Rayleigh-Taylor instability, P/e, A = 0.75, IC: light >> heavy.)"; }
};
using prod_A075L = keyword< prod_A075L_info, pegtl_string_t("prod_A075L") >;

struct hydroproductions_info {
  static std::string name() { return "P/eps"; }
  static std::string shortDescription() { return
    R"(Set MixMassFractionBeta SDE parameter(s) productions)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of strings used to parameterize
    the system of stochastic differential equations, configured in a
    mixmassfracbeta ... end block. Within the mixmassfracbeta ... end block the
    coefficients policy must be set to 'hydrotimescale' in order for the
    hydroproductions ... end  block to be in effect. The 'hydroproductions'
    keyword is then used to specify a list of strings, each specifying which
    turbulent kinetic energy production dividied by the dissipation rate (P/eps)
    data (from direct numerical simulations) should be used for the particular
    component integrated. Available P/eps data are defined in
    src/DiffEq/HydroProductions.h. Example: "productions prod_A05S prod_A05H
    prod_A05L prod_A05S prod_A05S end", which
    configures five P/eps data sets associated to 5 components, i.e., 5 scalar
    stochastic differential equations, integrated, specified and configured
    within the given mixmassfracbeta ... end block. The length of the
    hydroproductions vector depends on the number of scalar components and is
    controlled by the preceding keyword 'ncomp'. For mixmassfracbeta, ncomp is
    the actual number of scalar components * 4, since mixmassfractionbeta always
    computes 4 additional derived stochastic variables (in a diagnostic)
    fashion. See also MixMassFractionBeta::derived() in
    src/DiffEq/MixMassFractionBeta.h. Keywords allowed in a hydroproductions
    ... end block: )" + std::string("\'")
    + prod_A005H::string() + "\', \'"
    + prod_A005S::string() + "\', \'"
    + prod_A005L::string() + "\', \'"
    + prod_A05H::string() + "\', \'"
    + prod_A05S::string() + "\', \'"
    + prod_A05L::string() + "\', \'"
    + prod_A075H::string() + "\', \'"
    + prod_A075S::string() + "\', \'"
    + prod_A075L::string() + "\'. "
    + R"(For an example hydroproductions ... end block, see
      doc/html/walker_example_mixmassfracbeta.html.)"; }
  struct expect {
    using type = std::string;
    static std::string description() { return "string(s)"; }
  };
};
using hydroproductions =
  keyword< hydroproductions_info, pegtl_string_t("hydroproductions") >;

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
using gamma = keyword< gamma_info, pegtl_string_t("gamma") >;

struct ornstein_uhlenbeck_info {
  static std::string name() { return "Ornstein-Uhlenbeck"; }
  static std::string shortDescription() { return
    "Introduce the Ornstein-Uhlenbeck SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an ornstein-uhlenbeck ... end block,
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
  keyword< ornstein_uhlenbeck_info, pegtl_string_t("ornstein-uhlenbeck") >;

struct diag_ou_info {
  static std::string name() { return "Diagonal Ornstein-Uhlenbeck"; }
  static std::string shortDescription() { return
    "Introduce the diagonal Ornstein-Uhlenbeck SDE input block"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a diag_ou ... end
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
using diag_ou = keyword< diag_ou_info, pegtl_string_t("diag_ou") >;

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
using control = keyword< control_info, pegtl_string_t("control") >;

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
using smallcrush = keyword< smallcrush_info, pegtl_string_t("smallcrush") >;

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
using crush = keyword< crush_info, pegtl_string_t("crush") >;

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
using bigcrush = keyword< bigcrush_info, pegtl_string_t("bigcrush") >;

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
using verbose = keyword< verbose_info, pegtl_string_t("verbose") >;

struct benchmark_info {
  static std::string name() { return "benchmark"; }
  static std::string shortDescription() { return "Select benchmark mode"; }
  static std::string longDescription() { return
    R"(This keyword is used to select benchmark mode. In benchmark mode no large
       file output is performed, overriding the configuration in the control
       file.)";
  }
  using alias = Alias< b >;
};

using benchmark = keyword< benchmark_info, pegtl_string_t("benchmark") >;

struct feedback_info {
  static std::string name() { return "feedback"; }
  static std::string shortDescription() { return "Enable on-screen feedback"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable more detailed on-screen feedback on
       particular tasks and sub-tasks as they happen. This is useful for large
       problems and debugging.)";
  }
  using alias = Alias< f >;
};

using feedback = keyword< feedback_info, pegtl_string_t("feedback") >;

struct virtualization_info {
  static std::string name() { return "virtualization"; }
  static std::string shortDescription() { return
    R"(Set degree of virtualization)"; }
  static std::string longDescription() { return
    R"(This option is used to set the degree of virtualization
    (over-decomposition). The virtualization parameter is a real number
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
  keyword< virtualization_info, pegtl_string_t("virtualization") >;

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
using pdf = keyword< pdf_info, pegtl_string_t("pdf") >;

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
using stat = keyword< stat_info, pegtl_string_t("stat") >;

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
using input = keyword< input_info, pegtl_string_t("input") >;

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
using output = keyword< output_info, pegtl_string_t("output") >;

struct diagnostics_info {
  static std::string name() { return "diagnostics"; }
  static std::string shortDescription()
  { return "Specify the diagnostics file"; }
  static std::string longDescription() { return
    R"(This option is used to define the diagnostics file name.)";
  }
  using alias = Alias< d >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using diagnostics = keyword< diagnostics_info, pegtl_string_t("diagnostics") >;

struct reorder_info {
  static std::string name() { return "reorder"; }
  static std::string shortDescription() { return "Reorder mesh nodes"; }
  static std::string longDescription() { return
    R"(This option is used to instruct the mesh converter to not only convert
    but also reorder the mesh nodes using the advancing front technique.)";
  }
  using alias = Alias< r >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using reorder = keyword< reorder_info, pegtl_string_t("reorder") >;

struct group_info {
  static std::string name() { return "group"; }
  static std::string shortDescription() { return
    "Select test group(s) to run"; }
  static std::string longDescription() { return
    R"(This option can be used to select one or more test groups to run by
    specifying the full or a partial name of a test group. All tests of a
    selected group will be executed. If this option is not given, all test
    groups are executed by default. Examples: '--group make_list' - run only
    the 'make_list' test group, '--group Parser' - run the test groups that have
    the string 'Parser' in their name, e.g., groups 'Control/FileParser' and
    'Control/StringParser'.)";
  }
  using alias = Alias< g >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using group = keyword< group_info, pegtl_string_t("group") >;

struct inciter_info {
  static std::string name() { return "inciter"; }
  static std::string shortDescription() { return
    "Start configuration block for inciter"; }
  static std::string longDescription() { return
    R"(This keyword is used to select inciter. Inciter, is a continuum-realm
    shock hydrodynamics tool, solving a PDE.)";
  }
};
using inciter = keyword< inciter_info, pegtl_string_t("inciter") >;

struct user_defined_info {
  using code = Code< U >;
  static std::string name() { return "User-defined"; }
  static std::string shortDescription() { return
    "Select user-defined specification for a problem"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the user-define specification for a
    problem to be solved by a partial differential equation. The initial and
    boundary conditions are expected to be specified elsewhere in the input file
    to set up the problem. Example: "problem user_defined". This the default
    problem type.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};

using user_defined = keyword< user_defined_info, pegtl_string_t("user_defined") >;

struct shear_diff_info {
  using code = Code< S >;
  static std::string name() { return "Shear-diffusion"; }
  static std::string shortDescription() { return
    "Select the shear + diffusion test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the shear diffusion test problem. The
    initial and boundary conditions are specified to set up the test problem
    suitable to exercise and test the advection and diffusion terms of the
    scalar transport equation. Example: "problem shear_diff".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using shear_diff = keyword< shear_diff_info, pegtl_string_t("shear_diff") >;

struct dir_neu_info {
  using code = Code< D >;
  static std::string name() { return "Dirichlet & Neumann"; }
  static std::string shortDescription() { return
    "Select the Poisson equation test problem with Dirichlet and Neumann BCs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a Poisson equation test problem.
    Dirichlet and Neumann boundary conditions are simply hard-coded to set up
    the test problem, suitable to exercise and test the finite element
    discretization of the Laplace operator yielding both volume and boundary
    integral terms in its weak form. Example: "problem dir_neu".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using dir_neu = keyword< dir_neu_info, pegtl_string_t("dir_neu") >;

struct slot_cyl_info {
  using code = Code< Z >;
  static std::string name() { return "Zalesak's slotted cylinder"; }
  static std::string shortDescription() { return
    "Select Zalesak's slotted cylinder test problem"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Zalesak's slotted cylinder test
    problem. The initial and boundary conditions are specified to set up the
    test problem suitable to exercise and test the advection and diffusion
    terms of the scalar transport equation. Example: "problem slot_cyl".)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using slot_cyl = keyword< slot_cyl_info, pegtl_string_t("slot_cyl") >;

struct vortical_flow_info {
  using code = Code< V >;
  static std::string name() { return "Vortical flow"; }
  static std::string shortDescription() { return
    "Select the vortical flow test problem "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the vortical flow test problem. The
    purpose of this test problem is to test velocity errors generated by spatial
    operators in the presence of 3D vorticity and in particluar the
    superposition of planar and vortical flows, analogous to voritcity
    stretching. Example: "problem vortical_flow. For more details, see Waltz,
    et. al, "Manufactured solutions for the three-dimensional Euler equations
    with relevance to Inertial Confinement Fusion", Journal of Computational
    Physics 267 (2014) 196-209.)"; }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using vortical_flow =
  keyword< vortical_flow_info, pegtl_string_t("vortical_flow") >;

struct problem_info {
  using code = Code< r >;
  static std::string name() { return "problem"; }
  static std::string shortDescription() { return
    "Specify problem configuration for a partial differential equation solver";
  }
  static std::string longDescription() { return
    R"(This keyword is used to specify the problem configuration for a partial
    differential equation solver in the input file.)";
  }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + user_defined::string() + "\' | \'"
                  + shear_diff::string() + "\' | \'"
                  + dir_neu::string() + "\' | \'"
                  + vortical_flow::string() + "\' | \'"
                  + slot_cyl::string() + '\'';
    }
  };
};
using problem = keyword< problem_info, pegtl_string_t("problem") >;

struct compflow_navierstokes_info {
  using code = Code< N >;
  static std::string name() { return "Navier-Stokes"; }
  static std::string shortDescription() { return "Specify the Navier-Stokes "
    "(viscous) compressible flow physics configuration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Navier-Stokes (viscous) compressible "
    "flow physics configuration. Example: "compflow physics navierstokes end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using compflow_navierstokes =
  keyword< compflow_navierstokes_info, pegtl_string_t("navierstokes") >;

struct compflow_euler_info {
  using code = Code< E >;
  static std::string name() { return "Euler"; }
  static std::string shortDescription() { return "Specify the Euler (inviscid) "
    "compressible flow physics configuration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Euler (inviscid) compressible "
    "flow physics configuration. Example: "compflow physics euler end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using compflow_euler = keyword< compflow_euler_info, pegtl_string_t("euler") >;

struct advection_info {
  using code = Code< A >;
  static std::string name() { return "Advection"; }
  static std::string shortDescription() { return
    "Specify the advection physics configuration for a PDE "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the advection physics configuration for a
    PDE. Example: "transport physics advection end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using advection = keyword< advection_info, pegtl_string_t("advection") >;

struct advdiff_info {
  using code = Code< D >;
  static std::string name() { return "Advection + diffusion"; }
  static std::string shortDescription() { return
    "Specify the advection + diffusion physics configuration for a PDE "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the advection +diffusion physics
    configuration for a PDE. Example: "transport physics advdiff end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using advdiff = keyword< advdiff_info, pegtl_string_t("advdiff") >;

struct laplace_info {
  using code = Code< L >;
  static std::string name() { return "Laplace"; }
  static std::string shortDescription() { return
    "Specify the Laplace physics configuration for a PDE "; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Laplace physics configuration for a
    PDE. Example: "poisson physics laplace end")";
    }
  struct expect {
    static std::string description() { return "string"; }
  };
};
using laplace = keyword< laplace_info, pegtl_string_t("laplace") >;

struct physics_info {
  using code = Code< h >;
  static std::string name() { return "physics configuration"; }
  static std::string shortDescription() { return
    "Specify the physics configuration for a system of PDEs"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the physics configuration for a particular
    PDE system. Example: "physics navierstokes", which selects the Navier-Stokes
    equations for solving viscous compressible flow, given within the
    compflow ... end block. Valid options depend on the given block the keyword
    is used.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + advection::string() + "\' | \'"
                  + laplace::string() + "\' | \'"
                  + compflow_navierstokes::string() + "\' | \'"
                  + compflow_euler::string() + '\'';
    }
  };
};
using physics = keyword< physics_info, pegtl_string_t("physics") >;

struct pde_diffusivity_info {
  static std::string name() { return "diffusivity"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) diffusivity)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "diffusivity 5.0 2.0 3.0 end". The length of the vector depends on the
    particular type of PDE system and is controlled by the preceding keyword
    'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_diffusivity = keyword< pde_diffusivity_info, pegtl_string_t("diffusivity") >;

struct pde_lambda_info {
  static std::string name() { return "lambda"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) lambda)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "lambda 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of PDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_lambda = keyword< pde_lambda_info, pegtl_string_t("lambda") >;

struct pde_u0_info {
  static std::string name() { return "u0"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) u0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a vector of real numbers used to
    parameterize a system of partial differential equations. Example:
    "u0 5.0 2.0 3.0 end". The length of the vector depends on the particular
    type of PDE system and is controlled by the preceding keyword 'ncomp'.)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using pde_u0 = keyword< pde_u0_info, pegtl_string_t("u0") >;

struct pde_alpha_info {
  static std::string name() { return "alpha"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) alpha)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "alpha 5.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_alpha = keyword< pde_alpha_info, pegtl_string_t("alpha") >;

struct pde_beta_info {
  static std::string name() { return "beta"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) beta)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "beta 5.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_beta = keyword< pde_beta_info, pegtl_string_t("beta") >;

struct pde_p0_info {
  static std::string name() { return "p0"; }
  static std::string shortDescription() { return
    R"(Set PDE parameter(s) p0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify a real number used to
    parameterize a system of partial differential equations. Example:
    "p0 10.0".)"; }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real"; }
  };
};
using pde_p0 = keyword< pde_p0_info, pegtl_string_t("p0") >;

struct ctau_info {
  static std::string name() { return "ctau"; }
  static std::string shortDescription() { return
    R"(Set FCT mass diffusion coefficient, ctau)"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the mass diffusion coefficient used in
    flux-corrected transport, used for integrating transport equations. Example:
    "ctau 1.0".)"; }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static constexpr type upper = 1.0;
    static std::string description() { return "real"; }
    static std::string choices() {
      return "real between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "]";
    }
  };
};
using ctau = keyword< ctau_info, pegtl_string_t("ctau") >;

struct transport_info {
  static std::string name() { return "Transport"; }
  static std::string shortDescription() { return
    "Start configuration block for an transport equation"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an transport ... end block, used to
    specify the configuration for a transport equation type. Keywords allowed
    in an transport ... end block: )" + std::string("\'")
    + problem::string() + "\', \'"
    + physics::string() + "\', \'"
    + pde_diffusivity::string() + "\', \'"
    + pde_lambda::string() + "\', \'"
    + pde_u0::string() + "\'. "
    + R"(For an example transport ... end block, see
      doc/html/inicter_example_transport.html.)";
  }
};
using transport = keyword< transport_info, pegtl_string_t("transport") >;

struct poisson_info {
  static std::string name() { return "Poisson"; }
  static std::string shortDescription() { return
    "Start configuration block for a Poisson equation"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an poisson ... end block, used to
    specify the configuration for a partial differential equation of
    Poisson type. Keywords allowed in an poisson ... end block: )"
    + std::string("\'")
    + problem::string() + "\'. "
    + R"(For an example poisson ... end block, see
      doc/html/inicter_example_poisson.html.)";
  }
};
using poisson = keyword< poisson_info, pegtl_string_t("poisson") >;

struct sideset_info {
  static std::string name() { return "sideset"; }
  static std::string shortDescription() { return
    "Specify configuration for setting BC on a side set";
  }
  static std::string longDescription() { return
    R"(This keyword is used to specify boundary conditions on a side set for a
    solving partial differential equation.)";
  }
  struct expect {
    using type = std::string;
    static std::string description() { return "strings"; }
  };
};
using sideset = keyword< sideset_info, pegtl_string_t("sideset") >;

struct bc_dirichlet_info {
  static std::string name() { return "Dirichlet boundary condition"; }
  static std::string shortDescription() { return
    "Start configuration block describing Dirichlet boundary conditions"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce an bc_dirichlet ... end block, used to
    specify the configuration for setting Dirichlet boundary conditions for a
    partial differential equation. Keywords allowed in an bc_dirichlet ... end
    block: )" + std::string("\'")
    + sideset::string() + "\'. "
    + R"(For an example bc_dirichlet ... end block, see
      doc/html/inicter_example_poisson.html.)";
  }
};
using bc_dirichlet = keyword< bc_dirichlet_info, pegtl_string_t("bc_dirichlet") >;

struct id_info {
  static std::string name() { return "id"; }
  static std::string shortDescription() { return "ID"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify an ID, a positive integer.)";
  }
  struct expect {
    using type = uint64_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using id = keyword< id_info, pegtl_string_t("id") >;

struct mat_gamma_info {
  static std::string name() { return "gamma"; }
  static std::string shortDescription() { return "ratio of specific heats"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, ratio of specific
       heats.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_gamma = keyword< mat_gamma_info, pegtl_string_t("gamma") >;

struct mat_mu_info {
  static std::string name() { return "mu"; }
  static std::string shortDescription() { return "dynamic viscosity"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, dynamic
       viscosity.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_mu = keyword< mat_mu_info, pegtl_string_t("mu") >;

struct mat_cv_info {
  static std::string name() { return "cv"; }
  static std::string shortDescription() {
    return "specific heat at constant volume"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, specific heat at
       constant volume.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_cv = keyword< mat_cv_info, pegtl_string_t("cv") >;

struct mat_k_info {
  static std::string name() { return "k"; }
  static std::string shortDescription() { return "heat conductivity"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the material property, heat
       conductivity.)";
  }
  struct expect {
    using type = tk::real;
    static constexpr type lower = 0.0;
    static std::string description() { return "real"; }
  };
};
using mat_k = keyword< mat_k_info, pegtl_string_t("k") >;

struct material_info {
  static std::string name() { return "Material properties block"; }
  static std::string shortDescription() { return
    "Start configuration block for material properties"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a material ... end block, used to
    specify material properties. Keywords allowed in a material ... end
    block: )" + std::string("\'")
    + id::string()+ "\', \'"
    + mat_gamma::string()+ "\', \'"
    + mat_mu::string()+ "\', \'"
    + mat_cv::string()+ "\', \'"
    + mat_k::string() + "\'. "
    + R"(For an example material ... end block, see
      doc/html/inicter_example_compflow.html.)";
  }
};
using material = keyword< material_info, pegtl_string_t("material") >;

struct velocity_info {
  static std::string name() { return "velocity"; }
  static std::string shortDescription() { return
    "Specify velocity initial conditions";
  }
  static std::string longDescription() { return
    R"(This keyword is used to set initial conditions for the velocity field.)";
  }
  struct expect {
    using type = tk::real;
    static std::string description() { return "real(s)"; }
  };
};
using velocity = keyword< velocity_info, pegtl_string_t("velocity") >;

struct compflow_info {
  static std::string name() { return "Compressible flow"; }
  static std::string shortDescription() { return
    "Start configuration block for the compressible flow equations"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the compflow ... end block, used to
    specify the configuration for a system of partial differential equations,
    governing compressible fluid flow. Keywords allowed in an compflow ... end
    block: )" + std::string("\'")
    + problem::string() + "\'."
    + R"(For an example compflow ... end block, see
      doc/html/inicter_example_compflow.html.)";
  }
};
using compflow = keyword< compflow_info, pegtl_string_t("compflow") >;

struct rcb_info {
  static std::string name() { return "recursive coordinate bisection"; }
  static std::string shortDescription() { return
    "Select recursive coordinate bisection mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the recursive coordinate bisection (RCB)
    mesh partitioner. RCB is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.h for other valid options.)"; }
};
using rcb = keyword< rcb_info, pegtl_string_t("rcb") >;

struct rib_info {
  static std::string name() { return "recursive inertial bisection"; }
  static std::string shortDescription() { return
    "Select recursive inertial bisection mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the recursive inertial bisection (RIB)
    mesh partitioner. RIB is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.h for other valid options.)"; }
};
using rib = keyword< rib_info, pegtl_string_t("rib") >;

struct hsfc_info {
  static std::string name() { return "Hilbert space filling curve"; }
  static std::string shortDescription() { return
    "Select Hilbert Space Filling Curve (HSFC) mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the Hilbert Space Filling Curve (HSFC)
    mesh partitioner. HSFC is a geometry-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.h for other valid options.)"; }
};
using hsfc = keyword< hsfc_info, pegtl_string_t("hsfc") >;

struct mj_info {
  static std::string name() { return "multi-jagged"; }
  static std::string shortDescription() { return
    "Select multi-jagged (MJ) mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the multi-jagged (MJ) mesh partitioner.
    MJ is a geometry-based partitioner used to distribute an input mesh among
    processing elements. See
    Control/Options/PartitioningAlgorithm.h for other valid options.)"; }
};
using mj = keyword< mj_info, pegtl_string_t("mj") >;

struct phg_info {
  static std::string name() { return "hypergraph"; }
  static std::string shortDescription() { return
    "Select parallel hypergraph mesh partitioner"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the parallel hypergraph (PHG)
    mesh partitioner. PHG is a graph-based partitioner used to distribute an
    input mesh among processing elements. See
    Control/Options/PartitioningAlgorithm.h for other valid options.)"; }
};
using phg = keyword< phg_info, pegtl_string_t("phg") >;

struct algorithm_info {
  static std::string name() { return "algorithm"; }
  static std::string shortDescription() { return
    "Select mesh partitioning algorithm"; }
  static std::string longDescription() { return
    R"(This keyword is used to select a mesh partitioning algorithm. See
    Control/Options/PartitioningAlgorithm.h for valid options.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + rcb::string() + "\' | \'"
                  + rib::string() + "\' | \'"
                  + hsfc::string() + "\' | \'"
                  + mj::string() + "\' | \'"
                  + phg::string() + '\'';
    }
  };
};
using algorithm = keyword< algorithm_info, pegtl_string_t("algorithm") >;

struct partitioning_info {
  static std::string name() { return "partitioning"; }
  static std::string shortDescription() { return
    "Start configuration block for mesh partitioning"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce a partitioning ... end block, used to
    specify the configuration for mesh partitioning. Keywords allowed
    in a partitioning ... end block: )" + std::string("\'")
    + rcb::string() + "\'.";
  }
};
using partitioning = keyword< partitioning_info, pegtl_string_t("partitioning") >;

struct amr_uniform_info {
  static std::string name() { return "uniform"; }
  static std::string shortDescription() { return
    "Select uniform initial mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to select uniform initial mesh refinement.)"; }
};
using amr_uniform = keyword< amr_uniform_info, pegtl_string_t("uniform") >;

struct amr_initial_info {
  static std::string name() { return "initial refinement"; }
  static std::string shortDescription() { return
    "Configure initial mesh refinement (before t=0)"; }
  static std::string longDescription() { return
    R"(This keyword is used to select the type of initial mesh refinement that
    happens before t = 0. At this time, only uniform initial mehs refinement is
    supported.)"; }
  struct expect {
    static std::string description() { return "string"; }
    static std::string choices() {
      return '\'' + amr_uniform::string() + '\'';
    }
  };
};
using amr_initial = keyword< amr_initial_info, pegtl_string_t("initial") >;

struct amr_info {
  static std::string name() { return "AMR"; }
  static std::string shortDescription() { return
    "Start configuration block configuring adaptive mesh refinement"; }
  static std::string longDescription() { return
    R"(This keyword is used to introduce the amr ... end block, used to
    configure adaptive mesh refinement. Keywords allowed
    in this block: )" + std::string("\'")
    + amr_initial::string() + "\'.";
  }
};
using amr = keyword< amr_info, pegtl_string_t("amr") >;



////////// NOT YET FULLY DOCUMENTED //////////

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
using mix_iem = keyword<mix_iem_info,  pegtl_string_t("mix_iem") >;

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
using mix_iecm = keyword<mix_iecm_info,  pegtl_string_t("mix_iecm") >;

struct mix_dir_info {
  static std::string name() { return "Dirichlet"; }
  static std::string shortDescription() { return "Dirichlet"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_dir', is short for Dirichlet. It is a material
    mix model that explicitly satisfies the unit-sum requirement for all
    statistical samples.)";
  }
};
using mix_dir = keyword<mix_dir_info,  pegtl_string_t("mix_dir") >;

struct mix_gendir_info {
  static std::string name() { return "Generalized Dirichlet"; }
  static std::string shortDescription() { return "Generalized Dirichlet"; }
  static std::string longDescription() { return
    R"(Material mix model, 'mix_gendir', is short for Lochner's generalized
    Dirichlet. It is a material mix model that explicitly satisfies the
    unit-sum requirement for all statistical samples.)";
  }
};
using mix_gendir = keyword<mix_gendir_info,  pegtl_string_t("mix_gendir") >;

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
using hommix = keyword<hommix_info, pegtl_string_t("hommix") >;

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
using homhydro = keyword<homhydro_info,  pegtl_string_t("homhydro") >;

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
using homrt = keyword<homrt_info,  pegtl_string_t("homrt") >;

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
using spinsflow = keyword<spinsflow_info,  pegtl_string_t("spinsflow") >;

// This will go away once all the keywords below are documented
struct undefined_info {
  static std::string name() { return "undef"; }
  static std::string shortDescription() { return "undefined"; }
  static std::string longDescription() { return "Undefined."; }
};

using dmpi = keyword<undefined_info,  pegtl_string_t("dmpi") >;
using glbi = keyword<undefined_info,  pegtl_string_t("glbi") >;
using glob = keyword<undefined_info, pegtl_string_t("glob") >;
using pos_inviscid = keyword<undefined_info,  pegtl_string_t("pos_inviscid") >;
using pos_viscous = keyword<undefined_info,  pegtl_string_t("pos_viscous") >;
using mass_beta = keyword<undefined_info,  pegtl_string_t("mass_beta") >;
using hydro_slm = keyword<undefined_info,  pegtl_string_t("hydro_slm") >;
using hydro_glm = keyword<undefined_info,  pegtl_string_t("hydro_glm") >;
using mixrate_gamma = keyword<undefined_info,  pegtl_string_t("mixrate_gamma") >;
using freq_gamma = keyword<undefined_info,  pegtl_string_t("freq_gamma") >;
using position = keyword<undefined_info,  pegtl_string_t("position") >;
using hydro = keyword<undefined_info,  pegtl_string_t("hydro") >;
using mix = keyword<undefined_info,  pegtl_string_t("mix") >;
using nposition = keyword<undefined_info,  pegtl_string_t("nposition") >;
using ndensity = keyword<undefined_info,  pegtl_string_t("ndensity") >;
using nvelocity = keyword<undefined_info,  pegtl_string_t("nvelocity") >;
using nscalar = keyword<undefined_info,  pegtl_string_t("nscalar") >;
using nfreq = keyword<undefined_info,  pegtl_string_t("nfreq") >;
using SLM_C0 = keyword<undefined_info,  pegtl_string_t("C0") >;
using freq_gamma_C1 = keyword<undefined_info,  pegtl_string_t("C1") >;
using freq_gamma_C2 = keyword<undefined_info,  pegtl_string_t("C2") >;
using freq_gamma_C3 = keyword<undefined_info,  pegtl_string_t("C3") >;
using freq_gamma_C4 = keyword<undefined_info,  pegtl_string_t("C4") >;
using Beta_At = keyword<undefined_info,  pegtl_string_t("At") >;
using transported_scalar = keyword<undefined_info,  pegtl_string_t("Y") >;
using transported_scalar_fluctuation = keyword<undefined_info,  pegtl_string_t("y") >;
using velocity_x = keyword<undefined_info,  pegtl_string_t("U") >;
using velocity_fluctuation_x = keyword<undefined_info,  pegtl_string_t("u") >;
using velocity_y = keyword<undefined_info,  pegtl_string_t("V") >;
using velocity_fluctuation_y = keyword<undefined_info,  pegtl_string_t("v") >;
using velocity_z = keyword<undefined_info,  pegtl_string_t("W") >;
using velocity_fluctuation_z = keyword<undefined_info,  pegtl_string_t("w") >;
using pressure = keyword<undefined_info,  pegtl_string_t("P") >;
using pressure_fluctuation = keyword<undefined_info,  pegtl_string_t("p") >;
using density = keyword<undefined_info,  pegtl_string_t("R") >;
using density_fluctuation = keyword<undefined_info,  pegtl_string_t("r") >;

} // kw::

#endif // Keywords_h
