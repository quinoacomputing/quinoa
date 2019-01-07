################################################################################
#
# \file      cmake/FindCharm.cmake
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the Charm++
#
################################################################################

# Find Charm++
#
#  CHARM_FOUND        - True if the charmc compiler wrapper was found
#  CHARM_INCLUDE_DIRS - Charm++ include files paths
#
#  Set CHARM_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(CHARM_ROOT "/path/to/custom/charm") # prefer over system
#  find_package(Charm)
#  if(CHARM_FOUND)
#    # Link executables with the charmc wrapper
#    STRING(REGEX REPLACE "<CMAKE_CXX_COMPILER>" "${CHARM_COMPILER}"
#           CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE}")
#  endif()

function(_GET_CHARMINC _OUT_INC _charmc)
  file(STRINGS ${_charmc} _contents REGEX "^CHARMINC=")
  if(_contents)
    string(REGEX REPLACE "^CHARMINC=\"(.*)\"" "\\1" ${_OUT_INC} "${_contents}")
    set(${_OUT_INC} ${${_OUT_INC}} PARENT_SCOPE)
  else()
    message(FATAL_ERROR "file ${_charmc} does not exist")
  endif()
endfunction()

# If already in cache, be silent
if (CHARM_INCLUDE_DIRS AND CHARM_COMPILER AND CHARM_RUN)
  set (CHARM_FIND_QUIETLY TRUE)
endif()

INCLUDE(FindCygwin)

FIND_PROGRAM(CHARM_COMPILER
  NAMES charmc
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(CHARM_RUN
  NAMES charmrun
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_C_COMPILER
  NAMES ampicc
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_CXX_COMPILER
  NAMES ampicxx
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

FIND_PROGRAM(AMPI_RUN
  NAMES ampirun
  PATHS ${CHARM_ROOT}
        $ENV{CHARM_ROOT}
        ${CYGWIN_INSTALL_PATH}
        ${CMAKE_INSTALL_PREFIX}/charm
  PATH_SUFFIXES bin
)

if(CHARM_COMPILER)
  _GET_CHARMINC(HINTS_CHARMINC ${CHARM_COMPILER})
endif()

FIND_PATH(CHARM_INCLUDE_DIR NAMES charm.h
                            HINTS ${HINTS_CHARMINC}
                                  ${CHARM_ROOT}/include
                                  $ENV{CHARM_ROOT}/include
                                  ${CMAKE_INSTALL_PREFIX}/charm/include
                            PATH_SUFFIXES charm)
if(CHARM_INCLUDE_DIR)
  set(CHARM_INCLUDE_DIRS ${CHARM_INCLUDE_DIR})
else()
  set(CHARM_INCLUDE_DIRS "")
endif()

# Handle the QUIETLY and REQUIRED arguments and set CHARM_FOUND to TRUE if all
# listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Charm DEFAULT_MSG CHARM_COMPILER
                                  CHARM_INCLUDE_DIRS CHARM_RUN)

if(AMPI_C_COMPILER AND AMPI_CXX_COMPILER)
  set(ENABLE_AMPI true)
endif()

if(CHARM_COMPILER)
  include(CheckIncludeFiles)
  CHECK_INCLUDE_FILES("${CHARM_INCLUDE_DIR}/conv-mach-opt.h"
                      HAVE_CHARM_CONV_MACH_OPT)

  if (HAVE_CHARM_CONV_MACH_OPT)
    include(CheckSymbolExists)
    CHECK_SYMBOL_EXISTS(CMK_SMP "${CHARM_INCLUDE_DIR}/conv-mach-opt.h"
                        CHARM_SMP)
    if (CHARM_SMP)
      message(STATUS "Charm++ built in SMP mode")
    else()
      message(STATUS "Charm++ built in non-SMP mode")
    endif()
  endif()

endif()

MARK_AS_ADVANCED(CHARM_COMPILER CHARM_INCLUDE_DIRS CHARM_RUN ENABLE_AMPI
                 AMPI_C_COMPILER AMPI_CXX_COMPILER AMPI_RUN CHARM_SMP)
