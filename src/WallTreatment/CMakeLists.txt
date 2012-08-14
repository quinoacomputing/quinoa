cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
)

add_library(SDE
            Langevin.C
)

INSTALL(TARGETS SDE
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
