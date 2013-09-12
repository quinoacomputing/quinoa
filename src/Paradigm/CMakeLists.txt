cmake_minimum_required(VERSION 2.8.5)

project(Paradigm CXX)

include_directories(${QUINOA_SOURCE_DIR}/Paradigm
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Control
)

add_library(Paradigm
            Paradigm.C
)

INSTALL(TARGETS Paradigm
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
