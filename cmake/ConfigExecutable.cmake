################################################################################
#
# \file      ConfigExecutable.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2021 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Configure Charm++ executable targets
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

  # Add charmrun to install dir
  install(FILES ${CMAKE_BINARY_DIR}/charmrun
          PERMISSIONS OWNER_READ OWNER_EXECUTE OWNER_WRITE
                      GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE
          DESTINATION ${CMAKE_INSTALL_BINDIR})

  message(STATUS "Executable '${target}' configured")

endfunction()
