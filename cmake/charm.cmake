if(__addCharmModule)
  return()
endif()
set(__addCharmModule YES)

# Function 'addCharmModule' is used to add custom build commands and dependency
# for a Charm++ module
function(addCharmModule SRCPATH MODULE HOSTFILE PARTOF)
# Arguments:
#   SRCPATH:  Path where the .ci Charm++ module file resides, relative to
#             CMAKE_SOURCE_DIR.
#   MODULE:   Name of the Charm++ module.
#   HOSTFILE: Name of the .C file that includes the Charm++-generated .decl.h
#             and .def.h.
#   PARTOF:   Name of the library or executable the Charm++ module (and the
#             Charm++ host file) will link to. Note that the library or
#             executable PARTOF requires a custom link command as that must be
#             linked by the charmc wrapper.

  # Add custom command generating .decl.h and .def.h from .ci
  add_custom_command(OUTPUT ${MODULE}.decl.h ${MODULE}.def.h
                     DEPENDS ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci
                     COMMAND ${CHARM_COMPILER}
                             ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci)
  add_custom_target(${MODULE}CharmModule
                    DEPENDS ${MODULE}.decl.h ${MODULE}.def.h)
  add_dependencies(${PARTOF} ${MODULE}CharmModule)

  # Add (current) binary directory, i.e., which this function was called from,
  # to list of include directories as .decl.h and .def.h will be generated there
  include_directories(${PROJECT_BINARY_DIR})

  # Ignore warnings for HOSTFILE including the .decl.h and .def.h
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
  set_property(SOURCE ${HOSTFILE} APPEND_STRING PROPERTY COMPILE_FLAGS
               ${HOSTFILE_FLAGS})

  # Ignore warnings for linking of PARTOF
  set(PARTOF_FLAGS "-Wno-unused-parameter")
  set_property(TARGET ${PARTOF} APPEND_STRING PROPERTY LINK_FLAGS ${PARTOF_FLAGS})
  # Adding a RULE_LAUNCH_LINK prefixes the link line with a given command, which
  # is given here by a shell script that hijacks the link line, replaces the
  # linker by the charmc wrapper and links, see also the hijacking script.
  set_target_properties(${PARTOF} PROPERTIES RULE_LAUNCH_LINK
    "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/charmlink.sh ${CHARM_COMPILER}")

  # Echo status on adding Charm++ module
  message(STATUS "Add Charm++ module ${MODULE}.ci in ${HOSTFILE} at ${CMAKE_SOURCE_DIR}/${SRCPATH}/ linking to ${PARTOF}")
  message(STATUS "  Extra ${HOSTFILE} compiler flags: ${HOSTFILE_FLAGS}")
  message(STATUS "  Extra ${PARTOF} linker flags: ${PARTOF_FLAGS}")

endfunction(addCharmModule)
