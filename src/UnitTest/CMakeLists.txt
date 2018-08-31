cmake_minimum_required(VERSION 2.8.5)

project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            Assessment.C
            TUTSuite.C
            TUTTest.C
            MPIRunner.C)

target_include_directories(UnitTest PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Main
                           ${QUINOA_SOURCE_DIR}/UnitTest
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${TUT_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../LinSys
                           ${PROJECT_BINARY_DIR}/../LoadBalance
                           ${PROJECT_BINARY_DIR}/../UnitTest)

addCharmModule( "tuttest" "UnitTest" )
addCharmModule( "tutsuite" "UnitTest" )
addCharmModule( "mpirunner" "UnitTest" )
addCharmModule( "charmchild" "UnitTest" )
addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "testarray" "UnitTest" )
addCharmModule( "migrated_base" "UnitTest" )
addCharmModule( "migrated_inciter" "UnitTest" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.C).
add_dependencies("UnitTest" "unittestCharmModule")

# Add extra dependencies of UnitTest on those CharmModules that are required for
# testing Inciter.
if (ENABLE_INCITER)
  add_dependencies("UnitTest" "matcgCharmModule")
  add_dependencies("UnitTest" "diagcgCharmModule")
  add_dependencies("UnitTest" "distfctCharmModule")
  add_dependencies("UnitTest" "dgCharmModule")
  add_dependencies("UnitTest" "discretizationCharmModule")
  add_dependencies("UnitTest" "solverCharmModule")
  add_dependencies("UnitTest" "transporterCharmModule")
endif()

# Add extra dependency of TUTSuite on charm modules testing advanced array
# element placement using maps. This is required since TUTSuite spawns chare
# arrays testing advanced array element placement maps.
add_dependencies( "tutsuiteCharmModule" "linearmapCharmModule" )
add_dependencies( "tutsuiteCharmModule" "unsmeshmapCharmModule" )

set_target_properties(UnitTest PROPERTIES LIBRARY_OUTPUT_NAME quinoa_unittest)

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
