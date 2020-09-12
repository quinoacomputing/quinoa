################################################################################
#
# \file      charm.cmake
# \copyright 2012-2015 J. Bakosi,
#            2016-2018 Los Alamos National Security, LLC.,
#            2019-2020 Triad National Security, LLC.
#            All rights reserved. See the LICENSE file for details.
# \brief     Function used to setup a Charm++ module
#
################################################################################

if(__addCharmModule)
  return()
endif()
set(__addCharmModule YES)

# Function 'addCharmModule' is used to add custom build commands and dependency
# for a Charm++ module
function(addCharmModule MODULE PARTOF)
  # Arguments:
  #   MODULE:    Name of the Charm++ module.
  #   PARTOF:    Name of the library or executable the Charm++ module (and the
  #              Charm++ host file) will link to. Note that the library or
  #              executable PARTOF requires a custom link command as that must
  #              be linked by the charmc wrapper.

  # Add custom command generating .decl.h and .def.h from .ci
  add_custom_command(OUTPUT ${MODULE}.decl.h ${MODULE}.def.h
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.ci
                     COMMAND ${CHARM_COMPILER}
                             ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE}.ci)
  # Add custom target dependency for Charm++ module
  add_custom_target(${MODULE}CharmModule
                    DEPENDS ${MODULE}.decl.h ${MODULE}.def.h)
  # Add dependency of PARTOF on new Charm++ module
  add_dependencies(${PARTOF} ${MODULE}CharmModule)

endfunction(addCharmModule)
