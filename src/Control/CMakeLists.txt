cmake_minimum_required(VERSION 2.8.5)

project(QuinoaControl CXX)

include(charm)

if(HAS_MKL)
  set (MKLRNG_SOURCES "Options/MKLUniformMethod.C"
                      "Options/MKLGaussianMethod.C")
endif()

add_library(QuinoaControl
            FileParser.C
            StringParser.C
            Quinoa/InputDeck/Parser.C
            Quinoa/CmdLine/Parser.C
            Options/RNG.C
            ${MKLRNG_SOURCES}
)

INSTALL(TARGETS QuinoaControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)

project(RNGTestControl CXX)

add_library(RNGTestControl
            FileParser.C
            StringParser.C
            RNGTest/InputDeck/Parser.C
            RNGTest/CmdLine/Parser.C
            Options/RNG.C
            ${MKLRNG_SOURCES}
)

removeWarnings( "RNGTest/InputDeck/Parser.C;RNGTest/CmdLine/Parser.C;Options/RNG.C" )

INSTALL(TARGETS RNGTestControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)

project(MeshConvControl CXX)

add_library(MeshConvControl
            StringParser.C
            Options/RNG.C
            MeshConv/CmdLine/Parser.C
)

INSTALL(TARGETS MeshConvControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
