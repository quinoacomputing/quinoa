################################################################################
#
# \file      cmake/FindBoostMPLCartesianProduct.cmake
# \author    C. Junghans
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Find the BoostMPLCartesianProduct header
# \date      Mon 05 Dec 2016 14:21:53 AM MST
#
################################################################################

# Find BoostMPLCartesianProduct headers and libraries
#
#  BoostMPLCartesianProduct_FOUND        - True if BoostMPLCartesianProduct is found
#  BoostMPLCartesianProduct_INCLUDE_DIR  - BoostMPLCartesianProduct include files paths
#
#  Set BoostMPLCartesianProduct_ROOT before calling find_package to a path to add an additional
#  search path, e.g.,
#
#  Usage:
#
#  set(BoostMPLCartesianProduct_ROOT "/path/to/custom/boost/mpl") # prefer over system
#  find_package(BoostMPLCartesianProduct)
#  if(BoostMPLCartesianProduct_FOUND)
#    FindBoostMPLCartesianProduct.cmake(${BoostMPLCartesianProduct_INCLUDE_DIR})
#  endif()

# If already in cache, be silent
if (BoostMPLCartesianProduct_INCLUDE_DIR)
  set (BoostMPLCartesianProduct_FIND_QUIETLY TRUE)
endif()

# Look for the header file
FIND_PATH(BoostMPLCartesianProduct_INCLUDE_DIR NAMES boost/mpl/cartesian_product.hpp HINTS ${BoostMPLCartesianProduct_ROOT}/include)

# Handle the QUIETLY and REQUIRED arguments and set BoostMPLCartesianProduct_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(BoostMPLCartesianProduct DEFAULT_MSG BoostMPLCartesianProduct_INCLUDE_DIR)

MARK_AS_ADVANCED(BoostMPLCartesianProduct_INCLUDE_DIR)
