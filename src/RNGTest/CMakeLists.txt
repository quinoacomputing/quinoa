cmake_minimum_required(VERSION 2.8.5)

project(RNGTest CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(RNGTest
            TestU01.C
            TestU01Stack.C
            TestU01Suite.C
            SmallCrush.C
            Crush.C
            BigCrush.C
)

addCharmModule( "testu01" "RNGTest" )
addCharmModule( "testu01suite" "RNGTest" )
removeWarnings( "TestU01Suite.C;TestU01Stack.C;SmallCrush.C;Crush.C;BigCrush.C;TestU01.C" )

INSTALL(TARGETS RNGTest
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
