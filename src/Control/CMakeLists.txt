cmake_minimum_required(VERSION 2.8.5)

project(QuinoaControl CXX)

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

INSTALL(TARGETS RNGTestControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
