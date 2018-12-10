cmake_minimum_required(VERSION 2.8.5)

project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            Assessment.C
            TUTSuite.C
            TUTTest.C
            QuietCerr.C
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
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../LinSys
                           ${PROJECT_BINARY_DIR}/../LoadBalance
                           ${PROJECT_BINARY_DIR}/../UnitTest)

addCharmModule( "charmchild" "UnitTest" )
addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "testarray" "UnitTest" )
addCharmModule( "migrated_base" "UnitTest" )
addCharmModule( "quietcerr" "UnitTest" )

if (ENABLE_INCITER)
  addCharmModule( "migrated_inciter" "UnitTest" )
  addCharmModule( "tutsuiteinciter" "UnitTest" )
  addCharmModule( "tuttestinciter" "UnitTest" )
  addCharmModule( "mpirunnerinciter" "UnitTest" )
  add_dependencies( "tutsuiteinciterCharmModule" "mpirunnerinciterCharmModule" )
else()
  addCharmModule( "mpirunner" "UnitTest" )
  addCharmModule( "tutsuite" "UnitTest" )
  addCharmModule( "tuttest" "UnitTest" )
  add_dependencies( "tutsuiteCharmModule" "mpirunnerCharmModule" )
endif()

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.C).
if (ENABLE_INCITER)
  add_dependencies("UnitTest" "unittestinciterCharmModule")
else()
  add_dependencies("UnitTest" "unittestCharmModule")
endif()

# Add extra dependencies of UnitTest on those CharmModules that are required for
# testing Inciter.
if (ENABLE_INCITER)
  add_dependencies("UnitTest" "diagcgCharmModule"
                              "alecgCharmModule"
                              "distfctCharmModule"
                              "dgCharmModule"
                              "discretizationCharmModule"
                              "transporterCharmModule")
endif()

# Add extra dependency of TUTSuite on charm modules testing advanced array
# element placement using maps. This is required since TUTSuite spawns chare
# arrays testing advanced array element placement maps.
if (ENABLE_INCITER)
  add_dependencies( "tutsuiteinciterCharmModule"
                    "linearmapCharmModule" "unsmeshmapCharmModule" )
else()
  add_dependencies( "tutsuiteCharmModule"
                    "linearmapCharmModule" "unsmeshmapCharmModule" )
endif()

set_target_properties(UnitTest PROPERTIES LIBRARY_OUTPUT_NAME quinoa_unittest)

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
