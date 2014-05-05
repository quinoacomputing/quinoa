if(__addCharmModule)
	return()
endif()
set(__addCharmModule YES)


function(addCharmModule SRCPATH MODULE HOSTFILE PARTOF)
# Arguments:
#   SRCPATH:  path where the .ci Charm++ module file resides, relative to src,
#   i.e., CMAKE_SOURCE_DIR
#   MODULE:   name of the Charm++ module
#   HOSTFILE: name of the .C file (without the extension) that includes the
#             Charm++-generated .decl.h and .def.h
#   PARTOF:   name of the library or executable the Charm++ module will link to

  # Add custom command generating .decl.h and .def.h from .ci
  add_custom_command(OUTPUT ${MODULE}.decl.h ${MODULE}.def.h
                     DEPENDS ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci
                     COMMAND ${CHARM_COMPILER}
                             ${CMAKE_SOURCE_DIR}/${SRCPATH}/${MODULE}.ci)
  add_custom_target(${MODULE}CharmModule
                    DEPENDS ${MODULE}.decl.h ${MODULE}.def.h)
  add_dependencies(${PARTOF} ${MODULE}CharmModule)

  # Add (current) binary directory, i.e. which this function was called from, to
  # list of include directories as .decl.h and .def.h will be generated there
  include_directories(${PROJECT_BINARY_DIR})

  # Ignore some warnings for HOSTFILE including the Charm-generated .decl.h and
  # .def.h
  set_source_files_properties(${HOSTFILE}.C PROPERTIES COMPILE_FLAGS "-Wno-shadow -Wno-unused-variable -Wno-unused-parameter -Wno-deprecated-register -Wno-mismatched-tags -Wno-cast-align")

  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")  # intel-specific ignores
    set_source_files_properties(${HOSTFILE}.C PROPERTIES
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

  # Build list of compiler arguments with each include directory, -I..., -I...
  get_property(include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
               PROPERTY INCLUDE_DIRECTORIES)
  set(INCDIR_FLAGS "")
  foreach(dir ${include_dirs})
    message(STATUS "dir='${dir}'")
    LIST(APPEND INCDIR_FLAGS "-I${dir}")
  endforeach()

  # Add custom command compiling the Charm++ host file using the charmc wrapper
  add_custom_command(OUTPUT ${CMAKE_SOURCE_DIR}/${SRCPATH}/${HOSTFILE}.o
                     DEPENDS ${CMAKE_SOURCE_DIR}/${SRCPATH}/${HOSTFILE}.C
                     COMMAND ${CHARM_COMPILER}
                             ${CMAKE_CXX_FLAGS}
                             ${CXXFLAGS}
                             ${CMAKE_SHARED_MODULE_CXX_FLAGS}
                             ${INCDIR_FLAGS}
                             -o ${CMAKE_SOURCE_DIR}/${SRCPATH}/${HOSTFILE}.o
                             ${CMAKE_SOURCE_DIR}/${SRCPATH}/${HOSTFILE}.C)
  add_custom_target([CharmHost]${HOSTFILE}.o
                    DEPENDS ${CMAKE_SOURCE_DIR}/${SRCPATH}/${HOSTFILE}.o
                            ${MODULE}CharmModule)
  add_dependencies(${PARTOF} [CharmHost]${HOSTFILE}.o)

  message(STATUS "Add Charm++ module ${MODULE}.ci in ${HOSTFILE}.C at ${CMAKE_SOURCE_DIR}/${SRCPATH}/ linking to ${PARTOF}")

endfunction()
