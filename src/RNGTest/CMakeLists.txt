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

# Add extra dependency of RNGTest on rngtestCharmModule. This is required as one
# of the dependencies of RNGTest, e.g., TestU01Suite, refers to the main Charm++
# proxy defined in the Charm++ module rngtest (in Main/RNGTest.C).
add_dependencies("RNGTest" "rngtestCharmModule")

INSTALL(TARGETS RNGTest
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
