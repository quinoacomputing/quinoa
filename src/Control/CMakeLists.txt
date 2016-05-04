cmake_minimum_required(VERSION 2.8.5)

project(Control CXX)

#### Quinoa control build configuration ########################################

include(charm)

#### Inciter control ###########################################################
add_library(InciterControl
            FileParser.C
            StringParser.C
            Inciter/InputDeck/Parser.C
            Inciter/CmdLine/Parser.C
)

INSTALL(TARGETS InciterControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### RNGTest control ###########################################################
project(RNGTestControl CXX)

add_library(RNGTestControl
            FileParser.C
            StringParser.C
            RNGTest/InputDeck/Parser.C
            RNGTest/CmdLine/Parser.C
)


INSTALL(TARGETS RNGTestControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### UnitTest control ##########################################################
project(UnitTestControl CXX)

add_library(UnitTestControl
            StringParser.C
            FileParser.C
            UnitTest/CmdLine/Parser.C
)


INSTALL(TARGETS UnitTestControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### MeshConv control ##########################################################
project(MeshConvControl CXX)

add_library(MeshConvControl
            StringParser.C
            MeshConv/CmdLine/Parser.C
)

INSTALL(TARGETS MeshConvControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### Walker control ############################################################
add_library(WalkerControl
            FileParser.C
            StringParser.C
            Walker/InputDeck/Parser.C
            Walker/CmdLine/Parser.C
)

INSTALL(TARGETS WalkerControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

