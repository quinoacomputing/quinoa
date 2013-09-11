cmake_minimum_required(VERSION 2.8.5)

project(QuinoaControl CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/IO
                    ${TPL_DIR}/include/pegtl
)

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

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Control
                    ${TPL_DIR}/include
)

add_library(RNGTestControl
            Parser.C
            RNGTestParser.C
)

INSTALL(TARGETS RNGTestControl
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
