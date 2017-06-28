################################################################################
#
# \file      cmake/FindGmsh.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find Gmsh
#
################################################################################

# Find Gmsh
#
#  GMSH_FOUND - System has gmsh
#  GMSH_EXECUTABLE - The gmsh executable
#
#  Usage:
#
#  set(GMSH_ROOT "/path/to/custom/gmsh") # prefer over system
#  find_package(Gmsh)

if(GMSH_EXECUTABLE)
  # Already in cache, be silent
  set (GMSH_FIND_QUIETLY TRUE)
endif()

FIND_PATH(GMSH_EXECUTABLE NAMES gmsh)

# Handle the QUIETLY and REQUIRED arguments and set GMSH_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Gmsh DEFAULT_MSG GMSH_EXECUTABLE)

MARK_AS_ADVANCED(GMSH_EXECUTABLE)
