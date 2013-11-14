cmake_minimum_required(VERSION 2.8.5)

project(QuinoaControl CXX)

add_library(QuinoaControl
            FileParser.C
            StringParser.C
            Quinoa/InputDeck/Parser.C
            Quinoa/CmdLine/Parser.C
            Quinoa/Options/RNG.C
            Quinoa/Options/MKLUniformMethod.C
            Quinoa/Options/MKLGaussianMethod.C
            Quinoa/Options/Mass.C
            Quinoa/Options/Hydro.C
            Quinoa/Options/Mix.C
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
            Quinoa/Options/RNG.C
            Quinoa/Options/MKLUniformMethod.C
            Quinoa/Options/MKLGaussianMethod.C
)

INSTALL(TARGETS RNGTestControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
