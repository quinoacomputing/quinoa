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
  add_custom_target(${MODULE}CharmModule
                    DEPENDS ${MODULE}.decl.h ${MODULE}.def.h)
  add_dependencies(${PARTOF} ${MODULE}CharmModule)

  # Ignore warnings for linking of PARTOF
  set(PARTOF_FLAGS "-Wno-unused-parameter")
  set_property(TARGET ${PARTOF} APPEND_STRING PROPERTY LINK_FLAGS ${PARTOF_FLAGS})
  # Adding a RULE_LAUNCH_LINK prefixes the link line with a given command, which
  # is given here by a shell script that hijacks the link line, replaces the
  # linker by the charmc wrapper and links, see also the hijacking script.
  set_target_properties(${PARTOF} PROPERTIES RULE_LAUNCH_LINK
    "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/charmlink.sh ${CHARM_COMPILER}")

  # Echo status on adding Charm++ module
  message(STATUS "Add Charm++ module ${MODULE}.ci to ${CMAKE_CURRENT_SOURCE_DIR}/${PARTOF}")
  message(STATUS "  Extra linker flags: ${PARTOF_FLAGS}")

endfunction(addCharmModule)


function(removeWarnings HOSTFILES)
  # Arguments:
  #   HOSTFILES: List of .C files that directly or indirectly include the
  #              Charm++-generated .decl.h and .def.h.

  # Ignore warnings for HOSTFILES including the .decl.h and .def.h
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific ignores
    # Suppress the following Intel compiler warnings/remarks:
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
    # 3280: declaration hides member
    set(HOSTFILE_FLAGS "-diag-disable 177,181,1572,1599,1720,1944,2259,3280")
  else()
    set(HOSTFILE_FLAGS "-Wno-shadow -Wno-unused-variable -Wno-unused-parameter -Wno-deprecated-register -Wno-mismatched-tags -Wno-cast-align")
  endif()
  set_property(SOURCE ${HOSTFILES} APPEND_STRING PROPERTY COMPILE_FLAGS
               ${HOSTFILE_FLAGS})

  message(STATUS "  Host file(s): ${HOSTFILES}")
  message(STATUS "  Extra host file compiler flags: ${HOSTFILE_FLAGS}")

endfunction(removeWarnings)
