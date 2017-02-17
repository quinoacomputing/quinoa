cmake_minimum_required(VERSION 2.8.5)

project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            Assessment.C
            TUTSuite.C
            TUTTest.C
)

addCharmModule( "tuttest" "UnitTest" )
addCharmModule( "tutsuite" "UnitTest" )
addCharmModule( "charmchild" "UnitTest" )
addCharmModule( "migrated" "UnitTest" )
addCharmModule( "charmtimer" "UnitTest" )
addCharmModule( "testarray" "UnitTest" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.C).
add_dependencies("UnitTest" "unittestCharmModule")

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
