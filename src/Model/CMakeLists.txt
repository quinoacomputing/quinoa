cmake_minimum_required(VERSION 2.8.5)

project(Model CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(Model
            Model.C
            ModelException.C
)

INSTALL(TARGETS Model
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
