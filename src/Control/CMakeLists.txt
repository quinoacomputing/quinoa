cmake_minimum_required(VERSION 2.8.5)

project(QuinoaControl CXX)

add_library(QuinoaControl
            Parser.C
            QuinoaParser.C
)

INSTALL(TARGETS QuinoaControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)

project(RNGTestControl CXX)

add_library(RNGTestControl
            Parser.C
            RNGTestParser.C
)

INSTALL(TARGETS RNGTestControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
