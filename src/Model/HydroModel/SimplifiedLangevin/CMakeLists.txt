cmake_minimum_required(VERSION 2.8.5)

project(SimplifiedLangevin CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Model/HydroModel/SimplifiedLangevin
                    ${QUINOA_SOURCE_DIR}/Model/HydroModel
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(SimplifiedLangevin
            SimplifiedLangevin.C
)

INSTALL(TARGETS SimplifiedLangevin
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
