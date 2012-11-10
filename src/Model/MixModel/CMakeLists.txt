cmake_minimum_required(VERSION 2.8.5)

project(MixModel CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Model/MixModel
                    ${QUINOA_SOURCE_DIR}/Model
)

add_library(MixModel
            MixModel.C
)

INSTALL(TARGETS MixModel
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
