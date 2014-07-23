cmake_minimum_required(VERSION 2.8.5)

project(Control CXX)

#### Quinoa control build configuration ########################################

include(charm)

if(HAS_MKL)
  set (MKLRNG_SOURCES "Options/MKLUniformMethod.C"
                      "Options/MKLGaussianMethod.C")
endif()

#### RNGTest control ###########################################################
add_library(QuinoaControl
            FileParser.C
            StringParser.C
            Quinoa/InputDeck/Parser.C
            Quinoa/CmdLine/Parser.C
            Options/RNG.C
            ${MKLRNG_SOURCES}
)

removeWarnings("Options/RNG.C")

INSTALL(TARGETS QuinoaControl
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
            Options/RNG.C
            ${MKLRNG_SOURCES}
)

removeWarnings("RNGTest/InputDeck/Parser.C;RNGTest/CmdLine/Parser.C")

INSTALL(TARGETS RNGTestControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### UnitTest control ##########################################################
project(UnitTestControl CXX)

add_library(UnitTestControl
            StringParser.C
            UnitTest/CmdLine/Parser.C
)

removeWarnings("UnitTest/CmdLine/Parser.C")

INSTALL(TARGETS UnitTestControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

#### MeshConv control ##########################################################
project(MeshConvControl CXX)

add_library(MeshConvControl
            StringParser.C
            Options/RNG.C
            MeshConv/CmdLine/Parser.C
)

removeWarnings("MeshConv/CmdLine/Parser.C")

INSTALL(TARGETS MeshConvControl
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
