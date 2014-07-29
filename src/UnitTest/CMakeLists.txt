cmake_minimum_required(VERSION 2.8.5)

project(UnitTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(UnitTest
            TUTTest.C
            TUTSuite.C
)

addCharmModule( "tuttest" "UnitTest" )
addCharmModule( "tutsuite" "UnitTest" )
addCharmModule( "charmchild" "UnitTest" )

# Add extra dependency of UnitTest on unittestCharmModule. This is required as a
# dependency of UnitTest, TUTSuite, refers to the main Charm++ proxy defined in
# the Charm++ module unittest (in Main/UnitTest.C).
add_dependencies("UnitTest" "unittestCharmModule")

removeWarnings( "TUTSuite.C;TUTTest.C" )

INSTALL(TARGETS UnitTest
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
