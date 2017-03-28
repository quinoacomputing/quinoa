################################################################################
#
# \file      cmake/ConfigExecutable.cmake
# \author    J. Bakosi
# \copyright 2016, Los Alamos National Security, LLC.
# \brief     Configure Charm++ executable targets
# \date      Tue 28 Mar 2017 09:23:33 AM MDT
#
################################################################################


# ##############################################################################
# Function to configure executable targets
#
# config_executable( <target> )
#
# Mandatory arguments:
# --------------------
#
# target - Target on which to set properties
#
# Author: J. Bakosi
#
# ##############################################################################
function(config_executable target)

  # Add MPI compile flags
  if(MPI_COMPILE_FLAGS)
    set_target_properties(${target} PROPERTIES
                          COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
  endif()

  # Add MPI link flags
  if(MPI_LINK_FLAGS)
    set_target_properties(${target} PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
  endif()

  # Conditionally enforce static linking
  if(NOT BUILD_SHARED_LIBS)
    set_target_properties(${target} PROPERTIES LINK_FLAGS "-static")
  endif()

  INSTALL(TARGETS ${target}
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development)

  if(NOT CMAKE_GENERATOR STREQUAL "Ninja")
    add_custom_command(TARGET ${target} POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy
                               "${CMAKE_BINARY_DIR}/Main/charmrun"
                               "${CMAKE_BINARY_DIR}/charmrun")
  endif()

  message(STATUS "Executable '${target}' configured")

endfunction()
