if(__addCharmModule)
	return()
endif()
set(__addCharmModule YES)


function(addCharmModule SRCPATH MODULE HOSTFILE)
# Arguments:
#   SRCPATH:  path where the .ci Charm++ module file resides, relative to src/
#   MODULE:   name of the Charm++ module
#   HOSTFILE: name of the .C file that includes the Charm++-generated .decl.h
#             and .def.h

  # Add custom command generating .decl.h and .def.h from .ci
  add_custom_command(OUTPUT ${MODULE}.decl.h ${MODULE}.def.h
                     DEPENDS ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci
                     COMMAND ${CHARM_COMPILER}
                             ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci)

  add_custom_target(${MODULE}CharmModule
                    DEPENDS ${MODULE}.decl.h ${MODULE}.def.h)
  add_dependencies(${RNGTEST_EXECUTABLE} ${MODULE}CharmModule)

  # Ignore some warnings for the source file including the Charm-generated
  # .decl.h and .def.h
  set_source_files_properties(${HOSTFILE} PROPERTIES COMPILE_FLAGS "-Wno-shadow -Wno-unused-variable -Wno-unused-parameter -Wno-deprecated-register -Wno-mismatched-tags -Wno-cast-align")

  message(STATUS "Add Charm++ module ${MODULE}.ci in ${HOSTFILE} at ${CMAKE_SOURCE_DIR}/${SRCPATH}/")

endfunction()
