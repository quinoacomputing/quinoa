// *****************************************************************************
/*!
  \file      src/Control/Keywords.hpp
  \copyright 2012-2015 J. Bakosi,
             2016-2018 Los Alamos National Security, LLC.,
             2019-2021 Triad National Security, LLC.
             All rights reserved. See the LICENSE file for details.
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
    _kw::keyword_ template, defined in Control/Keyword.hpp. Specializing the
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
          static const type upper = 10;

          // Optional expected valid choices description, here giving
          // information on the expected type and the valid bounds. Note that
          // this can be any string, but if they exist, it is a good idea give
          // at least some information on the bounds, as below, since the bounds
          // are NOT displayed in the help for a keyword. This decision keeps
          // the bounds specifications generic since they can be any type. As a
          // result, the help structures, defined in HelpFactory.hpp, are simpler
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
  \see Control/Keyword.hpp
  \see Control/HelpFactory.hpp
*/
// *****************************************************************************
#ifndef Keywords_h
#define Keywords_h

#include <limits>

#undef I
#include <pegtl/contrib/alphabet.hpp>

#include "Types.hpp"
#include "Keyword.hpp"
#include "QuinoaBuildConfig.hpp"

//! Keywords used by all input deck and command line parsers
namespace kw {

using namespace tao::pegtl::alphabet;

struct end_info {
  static std::string name() { return "end"; }
  static std::string shortDescription() { return "End of an input block"; }
  static std::string longDescription() { return
    R"(The end of a block is given by the 'end' keyword in the input file.
    Example: "rngs ... end".)";
  }
};
using end = keyword< end_info, TAOCPP_PEGTL_STRING("end") >;

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
using help = keyword< help_info, TAOCPP_PEGTL_STRING("help") >;

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
using helpctr = keyword< helpctr_info, TAOCPP_PEGTL_STRING("helpctr") >;

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
using helpkw = keyword< helpkw_info, TAOCPP_PEGTL_STRING("helpkw") >;

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
using pdfs = keyword< pdfs_info, TAOCPP_PEGTL_STRING("pdfs") >;

struct raw_info {
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
    the init policies in DiffEq/InitPolicy.hpp for valid options.)"; }
};
using raw = keyword< raw_info, TAOCPP_PEGTL_STRING("raw") >;

struct zero_info {
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
    policies in DiffEq/InitPolicy.hpp for valid options.)"; }
};
using zero = keyword< zero_info, TAOCPP_PEGTL_STRING("zero") >;

struct ncomp_info {
  static std::string name() { return "ncomp"; }
  static std::string shortDescription() { return
    "Set number of scalar components for a system of differential equations"; }
  static std::string longDescription() { return
    R"(This keyword is used to specify the number of scalar
    components of a vector. 'ncomp' means "number of components". It is also
    used for specifying the number of scalar components of a transporter scalar
    (see also the keywords 'transport').)";
  }
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static std::string description() { return "uint"; }
  };
};
using ncomp = keyword< ncomp_info,  TAOCPP_PEGTL_STRING("ncomp") >;

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
using ttyi = keyword< ttyi_info, TAOCPP_PEGTL_STRING("ttyi") >;

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
using control = keyword< control_info, TAOCPP_PEGTL_STRING("control") >;

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
using verbose = keyword< verbose_info, TAOCPP_PEGTL_STRING("verbose") >;

struct charestate_info {
  static std::string name() { return "charestate"; }
  static std::string shortDescription() { return
    "Enable verbose chare state screen output"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable verbose Charm++ chare state collection and
    screen output. The chare state is displayed after a run is finished and the
    data collected is grouped by chare id (thisIndex), and within groups data
    is ordered by the time-stamp when a given chare member function is
    called. See src/Base/ChareState.hpp for details on what is collected. Note
    that to collect chare state, the given chare must be instrumented. Note that
    if quescence detection is enabled,
    chare state collection is also automatically enabled, but the chare state is
    only output if quiescence is detected (which also triggers an error).)";
  }
  using alias = Alias< S >;
};
using charestate = keyword< charestate_info, TAOCPP_PEGTL_STRING("state") >;

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

using benchmark = keyword< benchmark_info, TAOCPP_PEGTL_STRING("benchmark") >;

struct nonblocking_info {
  static std::string name() { return "nonblocking"; }
  static std::string shortDescription()
  { return "Select non-blocking migration"; }
  static std::string longDescription() { return
    R"(This keyword is used to select non-blocking, instead of the default
       blocking, migration. WARNING: This feature is experimental, not well
       tested, and may not always work as expected.)";
  }
  using alias = Alias< n >;
};

using nonblocking =
  keyword< nonblocking_info, TAOCPP_PEGTL_STRING("nonblocking") >;

struct lbfreq_info {
  static std::string name() { return "Load balancing frequency"; }
  static std::string shortDescription()
  { return "Set load-balancing frequency during time stepping"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the frequency of load-balancing during
       time stepping. The default is 1, which means that load balancing is
       initiated every time step. Note, however, that this does not necessarily
       mean that load balancing will be performed by the runtime system every
       time step, only that the Charm++ load-balancer is initiated. For more
       information, see the Charm++ manual.)";
  }
  using alias = Alias< l >;
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max()-1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using lbfreq = keyword< lbfreq_info, TAOCPP_PEGTL_STRING("lbfreq") >;

struct rsfreq_info {
  static std::string name() { return "Checkpoint/restart frequency"; }
  static std::string shortDescription()
  { return "Set checkpoint/restart frequency during time stepping"; }
  static std::string longDescription() { return
    R"(This keyword is used to set the frequency of dumping checkpoint/restart
       files during time stepping. The default is 1000, which means that
       checkpoint/restart files are dumped at every 1000th time step.)";
  }
  using alias = Alias< r >;
  struct expect {
    using type = std::size_t;
    static constexpr type lower = 1;
    static constexpr type upper = std::numeric_limits< type >::max()-1;
    static std::string description() { return "int"; }
    static std::string choices() {
      return "integer between [" + std::to_string(lower) + "..." +
             std::to_string(upper) + "] (both inclusive)";
    }
  };
};
using rsfreq = keyword< rsfreq_info, TAOCPP_PEGTL_STRING("rsfreq") >;

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
using feedback = keyword< feedback_info, TAOCPP_PEGTL_STRING("feedback") >;

struct version_info {
  static std::string name() { return "Show version"; }
  static std::string shortDescription() { return "Show version information"; }
  static std::string longDescription() { return
    R"(This keyword is used to display version information for the
       executable/tool on the standard output and exit successfully.)";
  }
  using alias = Alias< V >;
};
using version = keyword< version_info, TAOCPP_PEGTL_STRING("version") >;

struct license_info {
  static std::string name() { return "Show license"; }
  static std::string shortDescription() { return "Show license information"; }
  static std::string longDescription() { return
    R"(This keyword is used to display license information for the
       executable/tool on the standard output and exit successfully.)";
  }
  using alias = Alias< L >;
};
using license = keyword< license_info, TAOCPP_PEGTL_STRING("license") >;

struct trace_info {
  static std::string name() { return "trace"; }
  static std::string shortDescription()
  { return "Disable call and stack trace"; }
  static std::string longDescription() { return R"(This keyword can be used to
    disable the on-screen call trace and stack trace after an exception is
    thrown. Trace output is on by default and in some cases, the call and
    stack trace can be huge and not very helpful, hence this command line
    option.)"; }
  using alias = Alias< t >;
};
using trace = keyword< trace_info, TAOCPP_PEGTL_STRING("trace") >;

struct quiescence_info {
  static std::string name() { return "quiescence"; }
  static std::string shortDescription()
  { return "Enable quiescence detection"; }
  static std::string longDescription() { return
    R"(This keyword is used to enable the quiescence detection feature of
       Charm++, used to catch logic errors in the asynchronous control flow,
       resulting in deadlocks. This is useful for automated testing and
       debugging and does have some overhead, so it is off by default.)";
  }
  using alias = Alias< q >;
};
using quiescence =
  keyword< quiescence_info, TAOCPP_PEGTL_STRING("quiescence") >;

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
  keyword< virtualization_info, TAOCPP_PEGTL_STRING("virtualization") >;

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
using pdf = keyword< pdf_info, TAOCPP_PEGTL_STRING("pdf") >;

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
using stat = keyword< stat_info, TAOCPP_PEGTL_STRING("stat") >;

struct input_info {
  static std::string name() { return "input"; }
  static std::string shortDescription() { return "Specify the input file"; }
  static std::string longDescription() { return
    R"(This option is used to define the name of input file.)";
  }
  using alias = Alias< i >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using input = keyword< input_info, TAOCPP_PEGTL_STRING("input") >;

struct output_info {
  static std::string name() { return "output"; }
  static std::string shortDescription() { return "Specify the output file"; }
  static std::string longDescription() { return
    R"(This option is used to define the output file name. In MeshConv, this is
    used to specify the output mesh file name. In Inciter this is used to
    specify the output base filename. The base filename is appended by
    ".e-s.<meshid>.<numchares>.<chareid>", where 'e-s' probably stands for
    ExodusII sequence (the output file format), <meshid> counts the number of
    new meshes (this is incremented whenever the mesh is new compared to the
    previous iteration, due to, e.g., mesh refinement), <numchares> is the total
    number of mesh partitions, and <chareid> is the work unit (or mesh
    partition) id.)";
  }
  using alias = Alias< o >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using output = keyword< output_info, TAOCPP_PEGTL_STRING("output") >;

struct refined_info {
  static std::string name() { return "Refined field output"; }
  static std::string shortDescription() { return
    "Turn refined field output on/off"; }
  static std::string longDescription() { return
    R"(This keyword can be used to turn on/off refined field output, which
    refines the mesh and evaluates the solution on the refined mesh for saving
    the solution.)"; }
  struct expect {
    using type = bool;
    static std::string description() { return "string"; }
    static std::string choices() { return "true | false"; }
  };
};
using refined =keyword< refined_info, TAOCPP_PEGTL_STRING("refined") >;

struct screen_info {
  static std::string name() { return "screen"; }
  static std::string shortDescription() {
    return "Specify the screen output file"; }
  static std::string longDescription() { return
    R"(This option is used to set the screen output file name. The default is
    "<executable>_screen.log".)";
  }
  using alias = Alias< O >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using screen = keyword< screen_info, TAOCPP_PEGTL_STRING("screen") >;

struct restart_info {
  static std::string name() { return "checkpoint/restart directory name"; }
  static std::string shortDescription()
    { return "Specify the directory for restart files"; }
  static std::string longDescription() { return
    R"(This option is used to specify the directory name in which to save
    checkpoint/restart files.)";
  }
  using alias = Alias< R >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using restart = keyword< restart_info, TAOCPP_PEGTL_STRING("restart") >;

struct diagnostics_cmd_info {
  static std::string name() { return "diagnostics"; }
  static std::string shortDescription()
  { return "Specify the diagnostics file name"; }
  static std::string longDescription() { return
    R"(This option is used to define the diagnostics file name.)";
  }
  using alias = Alias< d >;
  struct expect {
    using type = std::string;
    static std::string description() { return "string"; }
  };
};
using diagnostics_cmd =
  keyword< diagnostics_cmd_info, TAOCPP_PEGTL_STRING("diagnostics") >;

struct reorder_cmd_info {
  static std::string name() { return "reorder"; }
  static std::string shortDescription() { return "Reorder mesh nodes"; }
  static std::string longDescription() { return
    R"(This keyword is used as a command line argument to instruct the mesh
    converter to not only convert but also reorder the mesh nodes using the
    advancing front technique. Reordering is optional in meshconv and
    inciter.)";
  }
  using alias = Alias< r >;
  struct expect {
    using type = bool;
    static std::string description() { return "string"; }
  };
};
using reorder_cmd = keyword< reorder_cmd_info, TAOCPP_PEGTL_STRING("reorder") >;

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
using group = keyword< group_info, TAOCPP_PEGTL_STRING("group") >;

struct inciter_info {
  static std::string name() { return "inciter"; }
  static std::string shortDescription() { return
    "Start configuration block for inciter"; }
  static std::string longDescription() { return
    R"(This keyword is used to select inciter. Inciter, is a continuum-realm
    shock hydrodynamics tool, solving a PDE.)";
  }
};
using inciter = keyword< inciter_info, TAOCPP_PEGTL_STRING("inciter") >;

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
using sideset = keyword< sideset_info, TAOCPP_PEGTL_STRING("sideset") >;

// This will go away once all the keywords below are documented
struct undefined_info {
  static std::string name() { return "undef"; }
  static std::string shortDescription() { return "undefined"; }
  static std::string longDescription() { return "Undefined."; }
};

} // kw::

#endif // Keywords_h
