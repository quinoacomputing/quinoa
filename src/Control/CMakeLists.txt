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

set_target_properties(InciterControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_incitercontrol)

INSTALL(TARGETS InciterControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

#### RNGTest control ###########################################################
project(RNGTestControl CXX)

add_library(RNGTestControl
            FileParser.C
            StringParser.C
            RNGTest/InputDeck/Parser.C
            RNGTest/CmdLine/Parser.C
)

set_target_properties(RNGTestControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_rngtestcontrol)

INSTALL(TARGETS RNGTestControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

#### UnitTest control ##########################################################
project(UnitTestControl CXX)

add_library(UnitTestControl
            StringParser.C
            FileParser.C
            UnitTest/CmdLine/Parser.C
)

set_target_properties(UnitTestControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_unittestcontrol)

INSTALL(TARGETS UnitTestControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

#### MeshConv control ##########################################################
project(MeshConvControl CXX)

add_library(MeshConvControl
            StringParser.C
            MeshConv/CmdLine/Parser.C
)

set_target_properties(MeshConvControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_meshconvcontrol)

INSTALL(TARGETS MeshConvControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

#### FileConv control ##########################################################
project(FileConvControl CXX)

add_library(FileConvControl
            StringParser.C
            FileConv/CmdLine/Parser.C
)

set_target_properties(FileConvControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_fileconvcontrol)

INSTALL(TARGETS FileConvControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

#### Walker control ############################################################
add_library(WalkerControl
            FileParser.C
            StringParser.C
            Walker/InputDeck/Parser.C
            Walker/CmdLine/Parser.C
)

set_target_properties(WalkerControl PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_walkercontrol)

INSTALL(TARGETS WalkerControl
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
