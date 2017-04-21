################################################################################
#
# \file      src/CMakeLists.txt
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Build quinoa
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
set(MINOR_VER 1)

# Set LA-CC number
set(LACC "16015")

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
  set(CMAKE_EXE_LINKER_FLAGS "-static")
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

# Set build type (will become compiler define in Main/QuinoaConfig.h)
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})
# Query regression tests direcetory (will become compiler define in
# Main/QuinoaConfig.h)
set(REGRESSION_DIR "${PROJECT_SOURCE_DIR}/../regression")

# Query target architecture
include(TargetArch)
target_architecture(ARCH)
message(STATUS "Target architecture: ${ARCH}")

# If not specified on the command-line, set default third-party libraries
# directory based on cmake-supplied compiler ID
if (NOT TPL_DIR)
  string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
  set(REL_TPL_DIR ${PROJECT_SOURCE_DIR}/../tpl/install/${COMP}-${ARCH})
  get_filename_component(ABS_TPL_DIR ${REL_TPL_DIR} ABSOLUTE)
  set(TPL_DIR ${ABS_TPL_DIR} CACHE STRING "Third-party libraries directory.")
  if (EXISTS ${TPL_DIR}/include)
    set(TPL_INCLUDE_DIR ${TPL_DIR}/include)
    message(STATUS "Third-party library dir: " ${TPL_DIR})
  endif()
endif()

# Find MPI and get its underlying compilers
include(MPICompilers)
get_mpi_compilers()

# Set compilers
set(COMPILER ${UNDERLYING_CXX_COMPILER})
set(MPI_COMPILER ${MPI_CXX_COMPILER})

# If we have met the mandatory MPI dependency, add a definition
if(${MPI_CXX_FOUND})
  add_definitions(-DQUINOA_CONFIG_MPI_ENABLED)
endif(${MPI_CXX_FOUND})

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

endif()

# Query build hostname
site_name(HOSTNAME)

# Grab current date to stick into the binary
execute_process(
  COMMAND "date" #"+%m/%d/%Y"
  OUTPUT_VARIABLE BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

# Set the requirement for the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
message(STATUS "Required language standard: C++${CMAKE_CXX_STANDARD}")

# Get git revision
include(GetGitRevisionDescription)
git_describe(GIT_SHA1 --tags)
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

  # Linker flags for clang
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)      # clear link flags
  set(CMAKE_EXE_LINKER_FLAGS "-Wno-missing-prototypes -Wno-unused-parameter")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Gnu will use libstdc++ by default. We do not offer libc++.

  # Compiler flags for gcc
  add_compiler_flag("-fdiagnostics-color")

  add_compiler_flag("-Wall")
  add_compiler_flag("-Wextra")
  add_compiler_flag("-pedantic")
  add_compiler_flag("-pedantic-errors")
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
  add_compiler_flag("-Wmissing-include-dirs")
  add_compiler_flag("-Wmissing-noreturn")
  add_compiler_flag("-Wpacked")
  add_compiler_flag("-Wpointer-arith")
  add_compiler_flag("-Wredundant-decls")
  add_compiler_flag("-Wshadow")
  add_compiler_flag("-Wstack-protector")
  add_compiler_flag("-Wstrict-aliasing=2")
  add_compiler_flag("-Wswitch-default")
  add_compiler_flag("-Wunreachable-code")
  add_compiler_flag("-Wunused")
  add_compiler_flag("-Wunused-parameter")
  add_compiler_flag("-Wvariadic-macros")
  add_compiler_flag("-Wwrite-strings")
  add_compiler_flag("-Wno-sign-compare")
  add_compiler_flag("-Wno-unused-function")
  add_compiler_flag("-Wno-stack-protector")

  # Linker flags for gcc
  set(CMAKE_EXE_LINKER_FLAGS "-Wno-unused-parameter")

  # Detect gcc version
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
                  OUTPUT_VARIABLE GCC_VERSION)
  string(REGEX MATCHALL "[0-9]+" GCC_VERSION_COMPONENTS ${GCC_VERSION})
  list(GET GCC_VERSION_COMPONENTS 0 GCC_MAJOR)
  list(GET GCC_VERSION_COMPONENTS 1 GCC_MINOR)
  message(STATUS "GCC version: ${GCC_MAJOR}.${GCC_MINOR}")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  # Include function that allows to set standard C++ library
  include(libstdcxx)

  # Offer switch between libc++ and libstdc++, using libc++ as default
  if (LIBCXX_FOUND)
    set_libstdcpp_vs_libcpp(CMAKE_CXX_FLAGS "libc++")
  endif()

  # Use the C++11 standard (CMAKE_CXX_STANDARD does not set this for intel)
  add_compiler_flag("-std=c++11")

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
  add_compiler_flag("-diag-disable 7,68,82,111,177,186,193,280,304,327,383,411,444,593,869,981,1418,1419,1572,2259,2407,2415,2547,2557,3280,3346,3656")

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
  set(CMAKE_EXE_LINKER_FLAGS "-diag-disable 869")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")  # pgi-specific settings

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")  # cray-specific settings

  # Use the C++11 standard (CMAKE_CXX_STANDARD does not set this for cray)
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
  # 111: indicated statement is not reachable
  # 186: unsigned integer is being compared to zero
  # 193: zero is used for the indicated undefined preprocessing identifier
  # 174: indicated expression has no effect
  # 177: variable was declared but never referenced
  # 181: argument is incompatible with the corresponding format string conversion
  # 261: access control is not specified ("private" by default)
  # 237: indicated selector expression is constant
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
  add_compiler_flag("-hnomessage=11709:261:12307:12489:193:826:737:997:3138:3137:7405:1418:177:181:401:11749:11548:1982:368:324:47:1:7:111:68:10153:186:550:21:237:174:7212")
  add_c_compiler_flag("-hnomessage=11709:301:3137:3138:7405:193:1977")

endif()

# Set runner configuration. This is used when generating code coverage and
# running regression tests using ctest. Note that ninja puts chamrun to <build>
# while make puts it into <build>/Main, but after any executable target we copy
# <build>/Main/charmrun to <build> for non-ninja builds. See also
# cmake/ConfigExecutable.cmake.
SET(RUNNER_REQUIRED false) # if true all executables must be prefixed by runner
if (DEFINED ENV{NERSC_HOST})
  SET(RUNNER "srun")
  SET(RUNNER_NCPUS_ARG "-n")
  SET(RUNNER_REQUIRED true)
elseif(ARCH MATCHES "ppc64")
  find_program(RUNNER ${CMAKE_CURRENT_SOURCE_DIR}/../script/cobalt-runjob.sh)
  SET(RUNNER_NCPUS_ARG "--np")
  SET(RUNNER_ARGS "")
  SET(RUNNER_REQUIRED true)
else()
  SET(RUNNER "${CMAKE_BINARY_DIR}/charmrun")
  SET(RUNNER_NCPUS_ARG "+p")
endif()

# Set default parallel/serial job runner extra arguments
if (NOT ARCH MATCHES "ppc64" AND
    NOT DEFINED RUNNER_ARGS AND
    NOT DEFINED ENV{CRAY_PRGENVPGI} AND
    NOT DEFINED ENV{CRAY_PRGENVGNU} AND
    NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
    NOT DEFINED ENV{CRAY_PRGENVINTEL})

  set(RUNNER_ARGS "--bind-to none --map-by node"
      CACHE STRING "Arguments to pass to the parallel/serial job runner, e.g., mpirun.")

endif()

# Echo runner configuration
message(STATUS "Parallel/serial job runner required: ${RUNNER_REQUIRED}")
message(STATUS "Parallel/serial job runner: ${RUNNER}")
message(STATUS "Parallel job runner ncpus arg: ${RUNNER_NCPUS_ARG}")
message(STATUS "Parallel/serial job runner args: ${RUNNER_ARGS}")

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

# Include cmake code to enable configuration for data layouts
include(ConfigureDataLayout)

# Include third-party libraries configuration
include(TPLs)

# Set names of executables
set(INCITER_EXECUTABLE inciter)
if (HAS_TESTU01)
  set(RNGTEST_EXECUTABLE rngtest)
  set(RNGTEST_SRC_DIR ${QUINOA_SOURCE_DIR}/RNGTest)
  set(RNGTEST_BIN_DIR ${PROJECT_BINARY_DIR}/RNGTest)
endif()
set(MESHCONV_EXECUTABLE meshconv)
set(WALKER_EXECUTABLE walker)
set(UNITTEST_EXECUTABLE unittest)

# Create set of all executables tested. These are the executables that do
# "interesting things", i.e., excluding unittest but keeping rngtest as that is
# (to be) regression tested
set(NONTEST_EXECUTABLES ${INCITER_EXECUTABLE}
                        ${RNGTEST_EXECUTABLE}
                        ${MESHCONV_EXECUTABLE}
                        ${WALKER_EXECUTABLE})

# Create set of all executables.
set(EXECUTABLES ${INCITER_EXECUTABLE}
                ${RNGTEST_EXECUTABLE}
                ${MESHCONV_EXECUTABLE}
                ${WALKER_EXECUTABLE}
                ${UNITTEST_EXECUTABLE})

# Make sure that the config file is in the search path
# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/QuinoaConfig.h.in"
  "${PROJECT_BINARY_DIR}/Main/QuinoaConfig.h"
)

# Remove build hostname from cmake cache
unset(HOSTNAME CACHE)

### Set RPATH so that shared libs paths stay in the binary after install ###
# Don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH false)

# When building, don't use the install RPATH
# (use later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH false)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
     "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../doc/quinoa.doxy
                 ${CMAKE_CURRENT_BINARY_DIR}/doc/quinoa.doxy @ONLY)
  add_custom_target(doc
                    ${DOXYGEN_EXECUTABLE}
                    ${CMAKE_CURRENT_BINARY_DIR}/doc/quinoa.doxy
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                    COMMENT "Documentation" VERBATIM)
  message(STATUS "Add target 'doc' to generate documentation to ./doc/html/index.html")
endif(DOXYGEN_FOUND)

include_directories(${QUINOA_SOURCE_DIR}
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/DiffEq
                    ${QUINOA_SOURCE_DIR}/PDE
                    ${QUINOA_SOURCE_DIR}/Walker
                    ${QUINOA_SOURCE_DIR}/Inciter
                    ${QUINOA_SOURCE_DIR}/LinSys
                    ${QUINOA_SOURCE_DIR}/LoadBalance
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Main
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/RNG
                    ${RNGTEST_SRC_DIR}
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/Particles
                    ${QUINOA_SOURCE_DIR}/UnitTest
                    ${MPI_CXX_INCLUDE_PATH}
                    ${MKL_INCLUDE_DIRS}
                    ${HDF5_INCLUDE_DIRS}
                    ${HYPRE_INCLUDE_DIRS}
                    ${PSTREAMS_INCLUDE_DIRS}
                    ${PUGIXML_INCLUDE_DIRS}
                    ${H5PART_INCLUDE_DIRS}
                    ${TUT_INCLUDE_DIRS}
                    ${CARTESIAN_PRODUCT_INCLUDE_DIRS}
                    ${PEGTL_INCLUDE_DIRS}
                    ${RANDOM123_INCLUDE_DIRS}
                    ${TESTU01_INCLUDE_DIRS}
                    ${RNGSSE2_INCLUDE_DIRS}
                    ${LAPACKE_INCLUDE_DIRS}
                    ${Zoltan2_INCLUDE_DIR}
                    ${SEACASExodus_INCLUDE_DIRS}
                    ${TPL_INCLUDE_DIR}
                    ${CHARM_INCLUDE_DIRS}
                    ${PROJECT_BINARY_DIR}/LoadBalance   # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Main          # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Walker        # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Inciter       # for Charm++ modules
                    ${RNGTEST_BIN_DIR}
                    ${PROJECT_BINARY_DIR}/UnitTest      # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/LinSys        # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/IO            # for Charm++ modules
)

# Components
add_subdirectory(Main)
add_subdirectory(UnitTest)
add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(DiffEq)
add_subdirectory(PDE)
add_subdirectory(Walker)
add_subdirectory(Inciter)
add_subdirectory(LoadBalance)
add_subdirectory(LinSys)
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(RNG)
if (HAS_TESTU01)
  add_subdirectory(RNGTest)
endif()
add_subdirectory(Statistics)
add_subdirectory(Particles)

# Conditionally enforce static link for all executables. This must be after all
# add_executable() calls creating the targets for executables in EXECUTABLES.
if (NOT BUILD_SHARED_LIBS)
  foreach(exe IN LISTS EXECUTABLES)
    set_property(TARGET ${exe} PROPERTY LINK_SEARCH_START_STATIC 1)
    set_property(TARGET ${exe} PROPERTY LINK_SEARCH_END_STATIC 1)
  endforeach()
endif()

# Optionally enable cppcheck static analysis
include(CppCheck)

# Enable cmake's ctest test suite used for regression testing
enable_testing()

# Include cmake code common to all regression tests
message(STATUS "Adding regression tests ...")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../regression
                 ${CMAKE_BINARY_DIR}/regression
                 EXCLUDE_FROM_ALL)

# Setup code coverage for unit tests
if(CODE_COVERAGE)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(unit doc/html unittest_coverage
                            ${RUNNER}
                            TESTRUNNER_ARGS
                            ${RUNNER_NCPUS_ARG} ${PROCESSOR_COUNT}
                            Main/${UNITTEST_EXECUTABLE} -v
                            DEPENDS ${NONTEST_EXECUTABLES}
                                    ${UNITTEST_EXECUTABLE})
endif()

# Setup code coverage for regression tests
if(CODE_COVERAGE)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(regression doc/html regression_coverage
                            ${CMAKE_CTEST_COMMAND}
                            TESTRUNNER_ARGS
                            -j${PROCESSOR_COUNT}
                            ${REGRESSION_COVERAGE_ARGS}
                            DEPENDS ${NONTEST_EXECUTABLES})
  if (REGRESSION_COVERAGE_ARGS)
    message(STATUS "Regression coverage test-runner arguments: '${REGRESSION_COVERAGE_ARGS}'")
  endif()
endif()

# Setup code coverage for all individual coverage targets included
if(CODE_COVERAGE)
  setup_target_for_all_coverage(all doc/html test_coverage
                                ${RUNNER} ${RUNNER_NCPUS_ARG}
                                DEPENDS ${NONTEST_EXECUTABLES}
                                        ${UNITTEST_EXECUTABLE})
endif()
