################################################################################
#
# \file      src/CMakeLists.txt
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019 Triad National Security, LLC.
#            All rights reserved See the LICENSE file for details.
# \brief     Build Quinoa
#
################################################################################

cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project(QUINOA C CXX)

set(CMAKE_VERBOSE_MAKEFILE 1)

# Quinoa cmake modules directory
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")

# Disallow in-source builds
include(DisallowInSourceBuilds)
disallow_in_source_builds()

# Set version number: major.minor
set(MAJOR_VER 0)
set(MINOR_VER 2)

# Set LA-CC number
set(LACC "17087")

# Sets CMAKE_INSTALL_{BIN,LIB}DIR needed for multilib support
include(GNUInstallDirs)

# Set CMAKE_BUILD_TYPE
include(BuildType)

# Set BUILD_SHARED_LIBS
include(BuildShared)

# Conditionally enforce static linking for executables added after this point.
# This must precede all add_executable() calls adding executables to link
# statically.
if (NOT BUILD_SHARED_LIBS)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
endif()

# Determine number of CPUs available
if (DEFINED ENV{SLURM_NPROCS})
  set(PROCESSOR_COUNT $ENV{SLURM_NPROCS})
else()
  include(ProcessorCount)
  ProcessorCount(PROCESSOR_COUNT)
  if(PROCESSOR_COUNT EQUAL 0)     # if ncpus could not be determined, assign 1
    set(PROCESSOR_COUNT "1")
  endif()
endif()
message(STATUS "Found ${PROCESSOR_COUNT} PEs")

# Define PE count for a long-running parallel regression tests
math(EXPR ManyPEs "${PROCESSOR_COUNT}")
if(PROCESSOR_COUNT GREATER 4)
  math(EXPR ManyPEs "${PROCESSOR_COUNT}/2")
endif()

# Convert CMAKE_BUILD_TYPE
string(TOUPPER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE)

# Set build type (will be exported to the source)
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})
# Query regression tests directory (will be exported to the source)
get_filename_component(REGRESSION_DIR
  "${PROJECT_SOURCE_DIR}/../tests/regression" ABSOLUTE)

# Query target architecture
include(TargetArch)
target_architecture(ARCH)
message(STATUS "Target architecture: ${ARCH}")

# Echo install dir
message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})

# Detect operating system type and version
include(DetectOS)
detect_os()
message(STATUS "Detected OS: ${HOST_OS}")
if(HOST_OS MATCHES "alpine")
  set(HOST_OS_ALPINE true)
endif()

# If not specified on the command-line, set default third-party libraries
# directory based on cmake-supplied compiler ID
if (NOT TPL_DIR)
  string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
  set(REL_TPL_DIR ${PROJECT_SOURCE_DIR}/../external/install/${COMP}-${ARCH})
  get_filename_component(ABS_TPL_DIR ${REL_TPL_DIR} ABSOLUTE)
  set(TPL_DIR ${ABS_TPL_DIR} CACHE STRING "Third-party libraries directory.")
endif()

if (EXISTS ${TPL_DIR}/include)
  set(TPL_INCLUDE_DIR ${TPL_DIR}/include)
endif()
message(STATUS "Third-party library dir: " ${TPL_DIR})

# Detect compiler version
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
    CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  include(DetectCompilerVersion)
endif()

set(ENABLE_TESTS true CACHE BOOL "Enable unit-, and regression tests, as well as building the unit test hardness.")
if (NOT ENABLE_TESTS)
  message(STATUS "Tests disabled.")
endif()

# Save contents of license file and copyright info in cmake variables
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE" LICENSE)
string(REGEX REPLACE "All rights reserved\.*" "All rights reserved. See --license for details." COPYRIGHT "${LICENSE}")
string(REGEX REPLACE ";" "\\\\n" COPYRIGHT "${COPYRIGHT}")
string(REGEX REPLACE ";" "\\\\n" LICENSE "${LICENSE}")
string(APPEND LICENSE "\\n\\nFor third-party licenses, see https://quinoacomputing.org/licenses.html")

# Include third-party libraries configuration
include(TPLs)

# Compilers
# =========
# There are two sets of compilers, given by the following CMake variables:
#
# 1. The underlying compiler suite
#    UNDERLYING_C_COMPILER
#    UNDERLYING_CXX_COMPILER
#
# 2. The MPI wrappers (of the underlying compilers)
#    MPI_C_COMPILER
#    MPI_CXX_COMPILER

# Find MPI and get its underlying compilers. If a traditional MPI (e.g.,
# OpenMPI) is used, the get_mpi_compilers cmake function uses the built-in CMake
# functionality to find the MPI wrappers and based on the wrappers detects the
# underyling compilers If Charm++'s Adaptive MPI (AMPI) is used, the Charm++
# TPL build builds Charm++ as well as AMPI and we use the AMPI wrappers to build
# MPI libraries.
if(NOT ENABLE_AMPI)
  include(MPICompilers)
  get_mpi_compilers()
  set(LINKER_COMPILER "${CHARM_COMPILER}")
else()
  set(MPI_C_COMPILER "${AMPI_C_COMPILER}")
  set(MPI_CXX_COMPILER "${AMPI_CXX_COMPILER}")
  set(UNDERLYING_C_COMPILER "${CMAKE_C_COMPILER}")
  set(UNDERLYING_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
  set(LINKER_COMPILER "${AMPI_CXX_COMPILER}")
  set(MPI_CXX_INCLUDE_DIRS "")
  MESSAGE(STATUS "AMPIrun: " ${AMPI_RUN})
endif()

MESSAGE(STATUS "MPI C compiler: " ${MPI_C_COMPILER})
MESSAGE(STATUS "MPI C++ compiler: " ${MPI_CXX_COMPILER})
MESSAGE(STATUS "Underlying C compiler: " ${UNDERLYING_C_COMPILER})
MESSAGE(STATUS "Underlying C++ compiler: " ${UNDERLYING_CXX_COMPILER})

MESSAGE(STATUS "Charmrun: " ${CHARM_RUN})
MESSAGE(STATUS "Linker compiler: " ${LINKER_COMPILER})

if (EXTRA_LINK_ARGS)
  MESSAGE(STATUS "Extra linker args (EXTRA_LINK_ARGS): ${EXTRA_LINK_ARGS}")
endif()

add_definitions(-DQUINOA_CONFIG_MPI_ENABLED)

option(ENABLE_AMR_TRACE OFF)

if(ENABLE_AMR_TRACE)
    add_definitions(-DENABLE_TRACE)
endif(ENABLE_AMR_TRACE)

# Set compilers
set(COMPILER ${UNDERLYING_CXX_COMPILER})
set(MPI_COMPILER ${MPI_CXX_COMPILER})

# Mac OS X specifics
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

  # On Mac OS X, set RPATH
  set(CMAKE_MACOSX_RPATH ON)

  # On Mac OS X, do not to resolve symbols when creating a shared library and
  # only resolve the symbols when linking executables. This behavior is
  # consistent with gcc, intel, and clang on linux.
  set(CMAKE_SHARED_LINKER_FLAGS
      "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")

  # On Mac OS X, due to changes in SnowLeopard (darwin10) linker changes, we add
  # "-no_compact_unwind" to the link line. This tells the linker to not produce
  # compact unwind information from the dwarf unwind info in .o files. Then at
  # runtime the darwin unwinder will fallback and use the slow dwarf unwind
  # info.
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} -Wl,-no_compact_unwind" )

  # Set minimum OSX SDK target version
  #set(CMAKE_OSX_DEPLOYMENT_TARGET 107 CACHE STRING "Minimum OS X version to target for deployment (at runtime); newer APIs weak linked. Set to empty string for default value." FORCE)

endif()

# Query build hostname
site_name(HOSTNAME)

# Grab current date to stick into the executables and doc (footer)
execute_process(COMMAND "date" OUTPUT_VARIABLE BUILD_DATE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

execute_process(COMMAND "date" "+%A, %b %d, %Y"
                OUTPUT_VARIABLE BUILD_DATE_SIMPLE)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE_SIMPLE "${BUILD_DATE_SIMPLE}")

# Set the requirement for the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Required language standard: C++${CMAKE_CXX_STANDARD}")

# Get git revision
include(GetGitRevisionDescription)
git_describe(GIT_SHA1 --tags)
git_rev_parse(GIT_SHA1_ONLY)
message(STATUS "git sha1: ${GIT_SHA1}")

# Macro to add C++ compiler flags
macro (add_compiler_flag flag)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
endmacro (add_compiler_flag)

# Macro to add C compiler flags
macro (add_c_compiler_flag flag)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}")
endmacro (add_c_compiler_flag)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")  # clang-specific settings

  # Include function that allows to set standard C++ library
  include(libstdcxx)

  # Offer switch between libc++ and libstdc++, using libc++ as default
  if (LIBCXX_FOUND)
    set_libstdcpp_vs_libcpp(CMAKE_CXX_FLAGS "libc++")
  endif()

  # Compiler flags for clang
  add_compiler_flag("-fcolor-diagnostics")

  add_compiler_flag("-Weverything")
  add_compiler_flag("-Wno-c++98-compat")
  add_compiler_flag("-Wno-c++98-compat-pedantic")
  add_compiler_flag("-Wno-padded")
  add_compiler_flag("-Wno-exit-time-destructors")
  add_compiler_flag("-Wno-global-constructors")
  add_compiler_flag("-Wno-weak-vtables")
  add_compiler_flag("-Wno-unknown-pragmas")
  add_compiler_flag("-Wno-unknown-warning-option")
  add_compiler_flag("-Wno-unused-command-line-argument")
  add_compiler_flag("-Wno-disabled-macro-expansion")
  add_compiler_flag("-Wno-reserved-id-macro")
  add_compiler_flag("-Wno-unknown-warning-option")

  # Linker flags for clang
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)      # clear link flags
  set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wno-missing-prototypes -Wno-unused-parameter")

  # Using clang on mac in release mode, use -O1, otherwise we get a runtime
  # error: auxlhs_complete() not found referenced by libquinoa_inciter.dylib
  if(CMAKE_SYSTEM_NAME MATCHES "Darwin" AND CMAKE_BUILD_TYPE MATCHES "REL")
    set(CMAKE_CXX_FLAGS_RELEASE "-O1 -DNDEBUG" CACHE STRING
        "CMAKE_CXX_FLAGS_RELEASE: Flags used by the compiler during release builds."
        FORCE)
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O1 -g -DNDEBUG" CACHE STRING
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO: Flags used by the compiler during release builds with debug info."
        FORCE)
  endif()

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Gnu will use libstdc++ by default. We do not offer libc++.

  # Compiler flags for gcc
  add_compiler_flag("-fdiagnostics-color")

  add_compiler_flag("-Wall")
  add_compiler_flag("-Wextra")
  add_compiler_flag("-Wcast-align")
  add_compiler_flag("-Wcast-qual")
  add_compiler_flag("-Wdisabled-optimization")
  add_compiler_flag("-Wfloat-equal")
  add_compiler_flag("-Wformat=2")
  add_compiler_flag("-Wformat-nonliteral")
  add_compiler_flag("-Wformat-security")
  add_compiler_flag("-Wformat-y2k")
  add_compiler_flag("-Wimport")
  add_compiler_flag("-Winit-self")
  add_compiler_flag("-Winvalid-pch")
  add_compiler_flag("-Wmissing-field-initializers")
  add_compiler_flag("-Wmissing-format-attribute")
  add_compiler_flag("-Wmissing-noreturn")
  add_compiler_flag("-Wpacked")
  add_compiler_flag("-Wpointer-arith")
  add_compiler_flag("-Wredundant-decls")
  add_compiler_flag("-Wshadow")
  add_compiler_flag("-Wstack-protector")
  add_compiler_flag("-Wstrict-aliasing=2")
  add_compiler_flag("-Wunreachable-code")
  add_compiler_flag("-Wunused")
  add_compiler_flag("-Wunused-parameter")
  add_compiler_flag("-Wvariadic-macros")
  add_compiler_flag("-Wwrite-strings")
  add_compiler_flag("-Wno-sign-compare")
  add_compiler_flag("-Wno-unused-function")
  add_compiler_flag("-Wno-stack-protector")
  add_compiler_flag("-Wno-expansion-to-defined")
  add_compiler_flag("-Wno-int-in-bool-context")
  add_compiler_flag("-Wno-cast-function-type")
  add_compiler_flag("-Wno-format-overflow")
  add_compiler_flag("-Wno-pragmas")

  # Linker flags for gcc
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wno-unused-parameter")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  # Include function that allows to set standard C++ library
  include(libstdcxx)

  # Offer switch between libc++ and libstdc++, using libstdc++ as default
  if (LIBCXX_FOUND)
    set_libstdcpp_vs_libcpp(CMAKE_CXX_FLAGS "libstdc++")
  endif()

  # Compiler flags for intel
  add_compiler_flag("-w3")       # enable diagnostics: remarks, warnings, errors
  add_compiler_flag("-Wremarks") # display remarks and comments

  # Globally suppress the following Intel compiler remarks:
  #     7: unrecognized token
  #    68: integer conversion resulted in a change of sign
  #    82: storage class is not first
  #   111: statement is unreachable
  #   186: pointless comparison of unsigned integer with zero
  #   177: function was declared but never referenced
  #   193: zero used for undefined preprocessing identifier
  #   280: selector expression is constant
  #   304: access control not specified ("public" by default)
  #   383: value copied to temporary, reference to temporary used
  #   327: NULL reference is not allowed
  #   411: defines no constructor for initializer_list<>
  #   444: destructor for base class is not virtual
  #   593: parameter was set but never used
  #   869: parameter was never referenced
  #   981: operands are evaluated in unspecified order
  #  1011: missing return statement at end of non-void function (constexpr)
  #  1292: unknown attribute "fallthrough"
  #  1418: external function definition with no prior declaration
  #  1419: external declaration in primary source file
  #  1572: floating-point equality and inequality comparisons are unreliable
  #  2259: non-pointer conversion from "int" to "char" may lose significant bits
  #  2407: the initialization of class will be done before that of member
  #  2415: variable of static storage duration was declared but never referenced
  #  2547: non-system include dir will be ignored when both a system and non-system include directory were specified
  #  2557: comparison between signed and unsigned operands
  #  3280: declaration hides member
  #  3346: dynamic exception specifications are deprecated
  #  3656: variable may be used before its value is set
  #  10120: overriding '-O3' with '-O2'
  add_compiler_flag("-diag-disable 7,68,82,111,177,186,193,280,304,327,383,411,444,593,869,981,1011,1292,1418,1419,1572,2259,2407,2415,2547,2557,3280,3346,3656,10120")

  # Compiler flags for optimization report
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=0 -qopt-report-phase=loop,vec -qopt-report-file=stderr")

  # No limit on inlining. This can increase the build time.
  if(CMAKE_BUILD_TYPE MATCHES "REL") # Release | RelWithDebInfo | MinSizeRel
    #add_compiler_flag("-no-inline-max-size")
    #add_compiler_flag("-no-inline-max-total-size")
    #add_compiler_flag("-no-inline-max-per-routine")
    add_compiler_flag("-diag-disable 11074,11076")
  endif()

  # Linker flags for intel
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -diag-disable 869")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")  # pgi-specific settings

  # Display all error messages (inform, warn, severe and fatal)
  add_compiler_flag("-Minform=inform")
  add_compiler_flag("-a")
  add_compiler_flag("--pedantic")
  add_compiler_flag("--display_error_number")

  # ignore warning on parameter was set but never used
  add_compiler_flag("--diag_suppress550")
  # ignore warning on variable was declared but never referenced
  add_compiler_flag("--diag_suppress177")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")  # cray-specific settings

  # Use the C++14 standard (CMAKE_CXX_STANDARD does not set this for cray)
  add_compiler_flag("-hstd=c++11")

  # Compiler flags for cray
  # enable errors, warnings, cautions, notes, comments
  add_compiler_flag("-hmsglevel_0")
  add_c_compiler_flag("-hmsglevel_0")

  # Globally suppress the following Cray compiler messages:
  # 1: source file does not end with a new-line character
  # 7: indicated token is not valid in this context
  # 21: type qualifiers are meaningless in this declaration
  # 47: macro has an incompatible redefinition
  # 68: integer conversion resulted in a change of sign
  # 82: storage class specifier is not the first specifier in the declaration
  # 111: indicated statement is not reachable
  # 186: unsigned integer is being compared to zero
  # 193: zero is used for the indicated undefined preprocessing identifier
  # 174: indicated expression has no effect
  # 177: variable was declared but never referenced
  # 181: argument is incompatible with the corresponding format string conversion
  # 261: access control is not specified ("private" by default)
  # 237: indicated selector expression is constant
  # 284: a "NULL" reference is not allowed
  # 301: indicated "typedef" name has already been declared (with same type)
  # 324: duplicate friend declaration
  # 368: class template defines no constructor to initialize ...
  # 401: destructor for base class is not virtual
  # 550: parameter "parts" is set but never used
  # 737: the using-declaration is ignored--it refers to the current namespace
  # 997: function is the virtual function override intended?
  # 826: parameter was never referenced
  # 1977: declaration hides variable
  # 1982: function ...::operator new has no corresponding member operator delete (to be called if an exception is thrown during initialization of an allocated object)
  # 3138: variable "this" assumed referenced by ...
  # 3137: variable "this" assumed modified by ...
  # 7212: variable ... is used before it is defined
  # 11709: diretory was specified as both a system and non-system include directory -- the non-system entry will be ignored
  # 12307: dynamic exception specifications are deprecated
  # 12489: constexpr non-static member function will not be implicitly 'const' in C++14
  # 7405: variable may be used before it is defined
  # 1418: a pointer expression is cast to a pointer with a greater alignment requirement
  # 10153: the initialization of class "Chare" will be done before that of member
  # 11749: declaration hides member
  # 11548: comparison between signed and unsigned operands
  add_compiler_flag("-hnomessage=11709:261:12307:12489:193:826:737:997:3138:3137:7405:1418:177:181:401:11749:11548:1982:368:324:47:1:7:111:68:10153:186:550:21:237:174:7212:284:82")
  add_c_compiler_flag("-hnomessage=11709:301:3137:3138:7405:193:1977")

else()

  message(WARNING "Compiler '${CMAKE_CXX_COMPILER_ID}' could not be identified.")

endif()

# Set runner configuration. This is used when generating code coverage and
# running regression tests using ctest. Note that ninja puts chamrun to <build>
# while make puts it into <build>/Main, but after any executable target we copy
# <build>/Main/charmrun to <build> for non-ninja builds. See also
# cmake/ConfigExecutable.cmake.

# Defaults for runner, used to run regression tests
set(RUNNER_REQUIRED false)
set(RUNNER "${CHARM_RUN}" CACHE STRING
    "Parallel/serial job runner, used to run regression tests")
set(RUNNER_NCPUS_ARG "+p" CACHE STRING
    "The argument used to specify the number of PEs for the runner")
set(RUNNER_ARGS "" CACHE STRING
    "Arguments to pass to the parallel/serial job runner, e.g., to mpirun.")
set(POSTFIX_RUNNER_ARGS "" CACHE STRING
    "Arguments to pass to the parallel/serial job runner at the end of the command line.")

# Host-, and/or architecture-specifc runner configuration
if (DEFINED ENV{NERSC_HOST})

  set(RUNNER "srun")
  set(RUNNER_NCPUS_ARG "-n")
  set(RUNNER_REQUIRED true)

elseif(DEFINED ENV{CRAY_PRGENVPGI} OR
       DEFINED ENV{CRAY_PRGENVGNU} OR
       DEFINED ENV{CRAY_PRGENVCRAY} OR
       DEFINED ENV{CRAY_PRGENVINTEL})

  set(RUNNER_ARGS "")

endif()

# Echo runner configuration
message(STATUS "Parallel/serial job runner required: ${RUNNER_REQUIRED}")
message(STATUS "Parallel/serial job runner: ${RUNNER}")
message(STATUS "Parallel job runner ncpus arg: ${RUNNER_NCPUS_ARG}")
message(STATUS "Parallel/serial job runner args: ${RUNNER_ARGS}")
message(STATUS "Parallel/serial job runner postfix args: ${POSTFIX_RUNNER_ARGS}")

# Find Python: required for code coverage (fastcov) and doc (m.css)
find_package(PythonInterp 3.6)

# Attempt to find tools required for code coverage analysis
include(DetectCodeCoverage)

# Set compiler flags to perform code coverage analysis
if(CODE_COVERAGE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
endif()

# Echo compiler flags
message(STATUS "C compiler flags: '${CMAKE_C_FLAGS}'")
message(STATUS "C++ compiler flags: '${CMAKE_CXX_FLAGS}'")
message(STATUS "Executable linker flags: '${CMAKE_EXE_LINKER_FLAGS}'")

# Include cmake code to enable configuration for data layouts
include(ConfigureDataLayout)

# Create set of all executables tested. These are the executables that do
# "interesting things", i.e., excluding unittest but keeping rngtest as that is
# (to be) regression tested
set(NONTEST_EXECUTABLES ${INCITER_EXECUTABLE}
                        ${RNGTEST_EXECUTABLE}
                        ${MESHCONV_EXECUTABLE}
                        ${WALKER_EXECUTABLE}
                        ${FILECONV_EXECUTABLE})

# Create set of all executables.
set(EXECUTABLES ${INCITER_EXECUTABLE}
                ${RNGTEST_EXECUTABLE}
                ${MESHCONV_EXECUTABLE}
                ${WALKER_EXECUTABLE}
                ${UNITTEST_EXECUTABLE}
                ${FILECONV_EXECUTABLE})


# Decide if the Exception class will write to std::cerr. This may not always be
# thread-safe, so it can be disabled here by setting this to false. The default
# is true, because it works most of the time.
if (NOT EXCEPTIONS_WRITE_TO_CERR)     # set default
  # will become compiler define in Main/QuinoaConfig.h
  set(EXCEPTIONS_WRITE_TO_CERR true CACHE BOOL "Exceptions will write std::cerr.")
endif()
message(STATUS "Exceptions write to cerr: ${EXCEPTIONS_WRITE_TO_CERR}")

# Configure cmake variables to pass to the build
configure_file( "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.hpp.in"
                "${PROJECT_BINARY_DIR}/Main/QuinoaConfig.hpp" )
configure_file( "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.cpp.in"
                "${PROJECT_BINARY_DIR}/Main/QuinoaConfig.cpp" ESCAPE_QUOTES)

# Remove build hostname from cmake cache
unset(HOSTNAME CACHE)

if (BUILD_SHARED_LIBS)

  ### Set RPATH so that shared libs paths stay in the binary after install ###
  # Don't skip the full RPATH for the build tree
  set(CMAKE_SKIP_BUILD_RPATH false)
  
  # When building, don't use the install RPATH
  # (use later on when installing)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH false)
  
  # Add the automatically determined parts of the RPATH
  # which point to directories outside the build tree to the install RPATH
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
  
  # the RPATH to be used when installing, but only if it's not a system directory
  list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
       "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
  if("${isSystemDir}" STREQUAL "-1")
     set(CMAKE_INSTALL_RPATH
       "${CMAKE_INSTALL_RPATH};${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
  endif("${isSystemDir}" STREQUAL "-1")

endif(BUILD_SHARED_LIBS)

# Components
if (CHARM_FOUND AND PSTREAMS_FOUND)
  add_subdirectory(Base)
else()
  PrintMissing(Base "CHARM_FOUND;PSTREAMS_FOUND")
endif()

if (CHARM_FOUND AND PEGTL_FOUND AND BRIGAND_FOUND)
  add_subdirectory(Control)
else()
  PrintMissing(Control "CHARM_FOUND;PEGTL_FOUND;BRIGAND_FOUND")
endif()

add_subdirectory(Main)

if (CHARM_FOUND AND PEGTL_FOUND AND BRIGAND_FOUND AND Zoltan2_FOUND)
  add_subdirectory(LoadBalance)
else()
  PrintMissing(LoadBalance "CHARM_FOUND;PEGTL_FOUND;BRIGAND_FOUND;Zoltan2_FOUND")
endif()

if (CHARM_FOUND AND PEGTL_FOUND AND BRIGAND_FOUND AND RANDOM123_FOUND)
  add_subdirectory(RNG)
else()
  PrintMissing(RNG "CHARM_FOUND;PEGTL_FOUND;BRIGAND_FOUND")
endif()

if (CHARM_FOUND AND PEGTL_FOUND AND BRIGAND_FOUND)
  add_subdirectory(Statistics)
else()
  PrintMissing(Statistics "CHARM_FOUND;PEGTL_FOUND;BRIGAND_FOUND")
endif()

if ((ENABLE_MESHCONV OR ENABLE_INCITER) AND BRIGAND_FOUND AND HIGHWAYHASH_FOUND)
  add_subdirectory(Mesh)
  add_subdirectory(IO)
endif()

if (ENABLE_INCITER AND HIGHWAYHASH_FOUND)
  add_subdirectory(PDE)
  add_subdirectory(Inciter)
endif()

if (ENABLE_UNITTEST AND ENABLE_TESTS)
  add_subdirectory(UnitTest)
endif()

if (ENABLE_WALKER)
  add_subdirectory(DiffEq)
  add_subdirectory(Walker)
endif()

if (ENABLE_RNGTEST)
  add_subdirectory(RNGTest)
endif()

# Conditionally enforce static link for all executables. This must be after all
# add_executable() calls creating the targets for executables in EXECUTABLES.
if (NOT BUILD_SHARED_LIBS)
  foreach(exe IN LISTS EXECUTABLES)
    set_property(TARGET ${exe} PROPERTY LINK_SEARCH_START_STATIC 1)
    set_property(TARGET ${exe} PROPERTY LINK_SEARCH_END_STATIC 1)
  endforeach()
endif()

# Add a target to generate documentation using doxygen and m.css
find_package(Doxygen 1.8.15)
find_package(MCSS)
if (DOXYGEN_FOUND AND MCSS_FOUND AND PYTHONINTERP_FOUND AND PYGMENTS_FOUND AND
    JINJA2_FOUND AND LATEX_FOUND)

  find_program( PERL perl )
  find_program( AWK awk )
  find_program( XARGS xargs )
  find_program( GIT git )
  find_program( HEAD head )
  find_program( RECENT_MODS
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/recent_mods.sh )
  find_program( EXTRACT_CMD_KEYWORDS
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/extract_cmd_keywords.sh )
  find_program( EXTRACT_CTR_KEYWORDS
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/extract_ctr_keywords.sh )

  if (PERL AND AWK AND XARGS AND GIT AND HEAD AND RECENT_MODS)

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../doc/quinoa.doxy
                   ${CMAKE_CURRENT_BINARY_DIR}/doc/quinoa.doxy @ONLY)

    file(GLOB doc_pages CONFIGURE_DEPENDS
         RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "../doc/pages/*.dox")
    foreach(page ${doc_pages})
      configure_file(${page} ${CMAKE_CURRENT_BINARY_DIR}/doc/${page} @ONLY)
    endforeach()

    # Generate targets <executable>_cmd_keywords for extracting and generating
    # reference documentation pages for command line keywords
    set(executables_cmd_keywords)
    foreach(exe IN LISTS EXECUTABLES)
      add_custom_target(${exe}_cmd_keywords
                        DEPENDS ${exe} ${doc_pages}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                        COMMAND ${EXTRACT_CMD_KEYWORDS} ${exe}
                        COMMENT "Extract ${exe}'s command-line keywords")
      list(APPEND executables_cmd_keywords "${exe}_cmd_keywords")
    endforeach()

    # Generate targets <executable>_ctr_keywords for extracting and generating
    # reference documentation pages for control file keywords
    set(executables_ctr_keywords)
    foreach(exe IN LISTS EXECUTABLES)
      add_custom_target(${exe}_ctr_keywords
                        DEPENDS ${exe} ${doc_pages}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                        COMMAND ${EXTRACT_CTR_KEYWORDS} ${exe}
                        COMMENT "Extract ${exe}'s control file keywords")
      list(APPEND executables_ctr_keywords "${exe}_ctr_keywords")
    endforeach()

    add_custom_command(OUTPUT xml
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc/quinoa.doxy
                               ${executables_cmd_keywords}
                               ${executables_ctr_keywords}
                       COMMAND ${RECENT_MODS} ${CMAKE_CURRENT_SOURCE_DIR}/..
                       COMMAND ${DOXYGEN_EXECUTABLE} quinoa.doxy
                       COMMAND ${MCSS_DOX2HTML5} --no-doxygen quinoa.doxy
                       VERBATIM USES_TERMINAL)

    add_custom_target(doc DEPENDS xml
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                      COMMENT "Documentation")

    message(STATUS "Add target 'doc' to generate documentation to ./doc/html/index.html")

  endif()

endif()

# Optionally enable cppcheck static analysis
include(CppCheck)

# Enable cmake's ctest test suite used for regression testing
enable_testing()

# Include cmake code common to all regression tests
if (ENABLE_TESTS)
  add_subdirectory(${REGRESSION_DIR}
                   ${CMAKE_BINARY_DIR}/tests
                   EXCLUDE_FROM_ALL)
endif()

# Setup code coverage for unit tests
if(CODE_COVERAGE AND ENABLE_TESTS)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(unit doc/html unittest_coverage ${RUNNER}
                            TESTRUNNER_ARGS ${RUNNER_NCPUS_ARG}
                                            ${PROCESSOR_COUNT}
                                            ${RUNNER_ARGS}
                                            Main/${UNITTEST_EXECUTABLE} -v
                            DEPENDS ${NONTEST_EXECUTABLES}
                                    ${UNITTEST_EXECUTABLE})
endif()

# Setup code coverage for regression tests
if(CODE_COVERAGE AND ENABLE_TESTS)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(regression doc/html regression_coverage
                            ${CMAKE_CTEST_COMMAND}
                            TESTRUNNER_ARGS -j${PROCESSOR_COUNT}
                                            ${REGRESSION_COVERAGE_ARGS} -LE insane
                            DEPENDS ${NONTEST_EXECUTABLES})
  if (REGRESSION_COVERAGE_ARGS)
    message(STATUS "Regression coverage test-runner arguments: '${REGRESSION_COVERAGE_ARGS}'")
  endif()
endif()

# Setup code coverage for all individual coverage targets included
if(CODE_COVERAGE AND ENABLE_TESTS)
  setup_target_for_all_coverage(all doc/html test_coverage
                                ${RUNNER} ${RUNNER_NCPUS_ARG}
                                DEPENDS ${NONTEST_EXECUTABLES}
                                        ${UNITTEST_EXECUTABLE})
endif()
