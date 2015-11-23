# Build quinoa

cmake_minimum_required(VERSION 2.8.5)

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

# Set release number
set(RELEASE "LA-CC-XX-XXX")

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "CMAKE_BUILD_TYPE not specified, setting to 'DEBUG'")
  set(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Build type. Possible values: DEBUG | RELEASE | RELWITHDEBINFO | MINSIZEREL" FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Echo build type
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

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

# Set build type
set(BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Build shared-libs by default
set(BUILD_SHARED_LIBS on)

# If not specified on the command-line, set default third-party libraries
# directory based on cmake-supplied compiler ID
if (NOT TPL_DIR)
  string(TOLOWER "${CMAKE_CXX_COMPILER_ID}" COMP)
  set(REL_TPL_DIR ${PROJECT_SOURCE_DIR}/../tpl/install/${COMP})
  get_filename_component(ABS_TPL_DIR ${REL_TPL_DIR} ABSOLUTE)
  set(TPL_DIR ${ABS_TPL_DIR} CACHE STRING "Third-party libraries directory.")
endif()

# Echo status on third-party libraries directory
if (TPL_DIR)
  MESSAGE(STATUS "Third-party library dir: " ${TPL_DIR})
else()
  message(FATAL_ERROR "Third-party library dir '${TPL_DIR}' is not specified")
endif()

# Add TPL_DIR/include to modules directory for TPLs that provide cmake
# FIND_PACKAGE code, such as Trilinos
SET(CMAKE_PREFIX_PATH ${TPL_DIR} ${CMAKE_PREFIX_PATH})

# Find MPI and get its underlying compilers
include(MPICompilers)
get_mpi_compilers()

# Set compilers
set(COMPILER ${UNDERLYING_CXX_COMPILER})
set(MPI_COMPILER ${MPI_CXX_COMPILER})
set(CHARM_COMPILER ${TPL_DIR}/charm/bin/charmc)
message(STATUS "Charm++ compiler wrapper (used to compile Charm++ interface (.ci) files and link Charm++ executables): " ${CHARM_COMPILER})

# Query build hostname
site_name(HOSTNAME)

# Grab current date to stick into the binary
execute_process(
  COMMAND "date" #"+%m/%d/%Y"
  OUTPUT_VARIABLE BUILD_DATE
)
string(REGEX REPLACE "[\r\n]" "" BUILD_DATE "${BUILD_DATE}")

# Set names of executables
set(INCITER_EXECUTABLE inciter)
set(RNGTEST_EXECUTABLE rngtest)
set(MESHCONV_EXECUTABLE meshconv)
set(WALKER_EXECUTABLE walker)
set(UNITTEST_EXECUTABLE unittest)

# Create set of all executables tested. These are the executables that do
# "interesting things", i.e., excluding unittest but keeping rngtest as that is
# (to be) regression tested
set(NONTEST_EXECUTABLES ${INCITER_EXECUTABLE} ${RNGTEST_EXECUTABLE}
                        ${MESHCONV_EXECUTABLE} ${WALKER_EXECUTABLE})

# Get git revision
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)

# Add compiler flags
include(CheckCXXCompilerFlag)

macro (add_compiler_flag flag)
  check_cxx_compiler_flag(${flag} COMPILER_HAS${flag})
  if(COMPILER_HAS${flag})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
  endif()
endmacro (add_compiler_flag)

# Try and use the following arguments for all compilers
add_compiler_flag("-Wall")
add_compiler_flag("-pedantic")
add_compiler_flag("-Wextra")
add_compiler_flag("-Wshadow")
add_compiler_flag("-Wformat")
add_compiler_flag("-Wsign-compare")
add_compiler_flag("-Wpointer-arith")
add_compiler_flag("-Winvalid-pch")
add_compiler_flag("-Wcast-align")
add_compiler_flag("-Wwrite-strings")
add_compiler_flag("-Wunused")
add_compiler_flag("-Woverloaded-virtual")
add_compiler_flag("-fdiagnostics-show-option")

#add_compiler_flag("-default-to-aout")

# Suppress warning for all compilers: ISO C++ forbids variable length array
add_compiler_flag("-Wno-vla")

# Require a C++11 compiler
check_cxx_compiler_flag(-std=c++11 COMPILER_HAS_CPP11)
check_cxx_compiler_flag(--c++11 COMPILER_HAS_CXX11)
if(COMPILER_HAS_CPP11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_HAS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --c++11")
else()
  MESSAGE(FATAL_ERROR "Need a compiler with C++11 support.")
endif()

# Include function that allows to set standard C++ library
include(libstdcxx)

# Compiler-specific settings
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")  # clang-specific settings

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fcolor-diagnostics")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Qunused-arguments")

  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

  # Experimental: use the C++14 standard
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

  # See all template error messages
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-backtrace-limit=0")

  # Don't warn on variable-length arrays
  add_compiler_flag("-Wno-vla-extension")
  # Don't warn on defined but unused functions
  add_compiler_flag("-Wno-unused-function")
  # Dont' warn on unused local typedef
  add_compiler_flag("-Wno-unused-local-typedef")
  # Warn on implicit conversions
  add_compiler_flag("-Wconversion")

  # Suppress the following warnings for Charm++ hostfiles including the
  # Charm-generated .decl.h and .def.h:
  set(HOSTFILE_FLAGS "-Wno-shadow -Wno-unused-variable -Wno-unused-parameter -Wno-deprecated-register -Wno-mismatched-tags -Wno-cast-align -Wno-zero-length-array -Wno-sign-compare -Wno-unused-private-field ")
  # Ignore warnings for linking with a Charm++ module
  set(PARTOF_FLAGS "-Wno-unused-parameter")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")  # gcc-specific settings

  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libstdc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

  # Experimental: use the C++14 standard
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

  add_compiler_flag("-Wno-literal-suffix")
  add_compiler_flag("-Wno-unused-local-typedefs")
  add_compiler_flag("-Wno-unused-function")
  add_compiler_flag("-Wno-deprecated-declarations")
  add_compiler_flag("-fdiagnostics-color")

  # Warn on implicit conversions
  add_compiler_flag("-Wconversion")
  # Dont' warn on multi-line comments
  add_compiler_flag("-Wno-comment")

  # Suppress the following warnings for Charm++ hostfiles including the
  # Charm-generated .decl.h and .def.h:
  set(HOSTFILE_FLAGS "-Wno-shadow -Wno-unused-variable -Wno-unused-parameter -Wno-cast-align -Wno-sign-compare -Wno-pedantic")
  # Ignore warnings for linking with a Charm++ module
  set(PARTOF_FLAGS "-Wno-unused-parameter")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific settings

  add_compiler_flag("-Wremarks")
  # Suppress the following Intel compiler remarks:
  #     7: unrecognized token
  #   111: statement is unreachable
  #   193: zero used for undefined preprocessing identifier
  #   304: access control not specified ("public" by default)
  #   367: duplicate friend declaration
  #   383: value copied to temporary, reference to temporary used
  #   411: defines no constructor for initializer_list<>
  #   444: destructor for base class is not virtual
  #   869: parameter was never referenced
  #   981: operands are evaluated in unspecified order
  #  1418: external function definition with no prior declaration
  #  2282: unrecognized GCC pragma
  #  2415: variable of static storage duration was declared but never referenced
  #  2547: path was specified as both a system and non-system include directory
  #  3346: dynamic exception specifications are deprecated
  #  3377: constexpr function return is non-constant
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 7,111,193,304,367,383,411,444,869,981,1418,2282,2415,2547,3346,3377 ")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=0 -qopt-report-phase=loop,vec -qopt-report-file=stderr")

  # No limit on inlining
  add_compiler_flag("-no-inline-max-size")
  add_compiler_flag("-no-inline-max-total-size")
  add_compiler_flag("-no-inline-max-per-routine")

  # Set libc++ as the default standard C++ library and allow selection by
  # default
  if (NOT STDLIBCPP)
    setlibstdcxx("libstdc++")
  else()
    setlibstdcxx(${STDLIBCPP})
  endif()

  # Suppress the following Intel compiler warnings/remarks for Charm++
  # hostfiles including the Charm-generated .decl.h and .def.h:
  #  177: variable was declared but never referenced
  #  181: argument of type "unsigned char *" is incompatible with format
  #       "%c", expecting argument of type "char *"
  # 1572: floating-point equality and inequality comparisons are unreliable
  # 1599: declaration hides variable
  # 1720: function "<class>::operator new" has no corresponding member
  #       operator delete (to be called if an exception is thrown during
  #       initialization of an allocated object)
  # 1944: declaration shadows a member of 'this'
  # 2259: non-pointer conversion from "int" to "unsigned char" may lose
  #       significant bits
  # 2557: comparison between signed and unsigned operands
  # 3280: declaration hides member
  set(HOSTFILE_FLAGS "-diag-disable 177,181,1572,1599,1720,1944,2259,2557,3280 ")
  # Ignore warnings for linking with a Charm++ module
  set(PARTOF_FLAGS "-Wno-unused-parameter")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")  # pgi-specific settings

endif()

# Echo Charm++ hostfile extra compiler flags
message(STATUS "Charm++ host extra compiler flags: ${HOSTFILE_FLAGS}")
# Echo Charm++ module extra link flags
message(STATUS "Charm++ module extra linker flags: ${PARTOF_FLAGS}")

# Set compiler flags for use of libc++. If nothing is given, do nothing; gnu
# and intel will use libstdc++, while clang will use libc++ by default.
if(LIBCXX_CPP)
   # Use libc++ (http://libcxx.llvm.org)
   # Use libc++abi (http://libcxxabi.llvm.org)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -lc++abi")
endif()

# Set default charmrun path. This is used when generating code coverage and
# running regression test.
SET(CHARMRUN "${CMAKE_BINARY_DIR}/Main/charmrun")
# Ninja puts chamrun to build dir
IF(CMAKE_GENERATOR STREQUAL "Ninja")
  SET(CHARMRUN "${CMAKE_BINARY_DIR}/charmrun")
ENDIF()
message(STATUS "charmrun executable: ${CHARMRUN}")

# Attempt to find tools required for code coverage analysis
include(DetectCodeCoverage)

# Set compiler flags to perform code coverage analysis
if(CODE_COVERAGE)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -O0 --coverage -fprofile-arcs -ftest-coverage")
endif()

# Set data layout for particle properties
# Available options
set(LAYOUT_VALUES "particle" "equation")
# Initialize all to off
set(LAYOUT_PARTICLE_MAJOR off)  # 0
set(LAYOUT_EQUATION_MAJOR off)  # 1
# Set default and select from list
set(LAYOUT "particle" CACHE STRING "Particle properties layout. Default: (particle-major). Available options: ${LAYOUT_VALUES}(-major).")
SET_PROPERTY (CACHE LAYOUT PROPERTY STRINGS ${LAYOUT_VALUES})
STRING (TOLOWER ${LAYOUT} LAYOUT)
LIST (FIND LAYOUT_VALUES ${LAYOUT} LAYOUT_INDEX)
# Evaluate selected option and put in a define for it
IF (${LAYOUT_INDEX} EQUAL 0)
  set(LAYOUT_PARTICLE_MAJOR on)
ELSEIF (${LAYOUT_INDEX} EQUAL 1)
  set(LAYOUT_EQUATION_MAJOR on)
ELSEIF (${LAYOUT_INDEX} EQUAL -1)
  MESSAGE(FATAL_ERROR "Particle properties layout '${LAYOUT}' not supported, valid entries are ${LAYOUT_VALUES}(-major).")
ENDIF()
message(STATUS "Particle properties layout (LAYOUT): " ${LAYOUT} "(-major)")

# Include third-party libraries configuration
include(TPLs)

# Make sure that the config file is in the search path
# Configure a header file to pass info to the code
configure_file(
  "${PROJECT_SOURCE_DIR}/Main/Config.h.in"
  "${PROJECT_BINARY_DIR}/Main/Config.h"
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
                    ${QUINOA_SOURCE_DIR}/Walker
                    ${QUINOA_SOURCE_DIR}/Inciter
                    ${QUINOA_SOURCE_DIR}/LinSys
                    ${QUINOA_SOURCE_DIR}/LoadBalance
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Main
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/RNG
                    ${QUINOA_SOURCE_DIR}/RNGTest
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/UnitTest
                    ${TPL_DIR}/include
                    ${TPL_DIR}/include/pegtl
                    ${TPL_DIR}/include/tut
                    ${TPL_DIR}/charm/include
                    ${PROJECT_BINARY_DIR}/LoadBalance   # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Main          # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Walker        # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/Inciter       # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/RNGTest       # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/UnitTest      # for Charm++ modules
                    ${PROJECT_BINARY_DIR}/LinSys        # for Charm++ modules
)

# Components
add_subdirectory(Main)
add_subdirectory(UnitTest)
add_subdirectory(Base)
add_subdirectory(Control)
add_subdirectory(DiffEq)
add_subdirectory(Walker)
add_subdirectory(Inciter)
add_subdirectory(LoadBalance)
add_subdirectory(LinSys)
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(RNG)
add_subdirectory(RNGTest)
add_subdirectory(Statistics)

# Enable cmake's ctest test suite used for regression testing
enable_testing()

# Add cmake code describing how to build numdiff
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../tool/src/ndiff
                 ${CMAKE_BINARY_DIR}/regression/ndiff)

# Set mpirun arguments controling process binding
set(MPIRUN_BIND_ARGS "-bind-to none -map-by node")
message(STATUS "mpirun process binding for ctest: ${MPIRUN_BIND_ARGS}")

# Include cmake code common to all regression tests
message(STATUS "Adding regression tests ...")
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../regression
                 ${CMAKE_BINARY_DIR}/regression
                 EXCLUDE_FROM_ALL)

# Setup code coverage for unit tests
if(CODE_COVERAGE)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(unit doc/html unittest_coverage
                            ${CHARMRUN}
                            TESTRUNNER_ARGS
                            +p${PROCESSOR_COUNT} Main/${UNITTEST_EXECUTABLE} -v
                            DEPENDS ${NONTEST_EXECUTABLES}
                                    ${UNITTEST_EXECUTABLE})
endif()

# Setup code coverage for regression tests
if(CODE_COVERAGE)
  # Setup test coverage target. Make it dependend on all quinoa executables.
  setup_target_for_coverage(regression doc/html regression_coverage
                            ${CMAKE_CTEST_COMMAND}
                            TESTRUNNER_ARGS -j${PROCESSOR_COUNT}
                            DEPENDS ${NONTEST_EXECUTABLES} numdiff)
endif()

# Setup code coverage for all individual coverage targets included
if(CODE_COVERAGE)
  setup_target_for_all_coverage(all doc/html test_coverage
                                DEPENDS ${NONTEST_EXECUTABLES}
                                        ${UNITTEST_EXECUTABLE}
                                        numdiff)
endif()
