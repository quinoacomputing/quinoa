cmake_minimum_required(VERSION 2.8.5)

project(RegTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(RegTest
            RegSuite.C
            Regression.C
)

addCharmModule( "regression" "RegTest" )
addCharmModule( "regsuite" "RegTest" )

# Add extra dependency of RegTest on regtestCharmModule. This is required as a
# dependency of RegTest, RegSuite, refers to the main Charm++ proxy defined in
# the Charm++ module regtest (in Main/RegTest.C).
add_dependencies("RegTest" "regtestCharmModule")

removeWarnings( "RegSuite.C;Regression.C" )

INSTALL(TARGETS RegTest
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
