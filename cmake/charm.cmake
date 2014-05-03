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

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific ignores
    set_source_files_properties(${HOSTFILE} PROPERTIES
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
      COMPILE_FLAGS "-diag-disable 177,181,1572,1599,1720,1944,2259,3280")
  endif()

  message(STATUS "Add Charm++ module ${MODULE}.ci in ${HOSTFILE} at ${CMAKE_SOURCE_DIR}/${SRCPATH}/")

endfunction()
