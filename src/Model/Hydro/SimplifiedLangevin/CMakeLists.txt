cmake_minimum_required(VERSION 2.8.5)

project(SimplifiedLangevin CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model/Hydro/SimplifiedLangevin
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Control
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
