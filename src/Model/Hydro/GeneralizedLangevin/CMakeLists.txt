cmake_minimum_required(VERSION 2.8.5)

project(GeneralizedLangevin CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model/Hydro/GeneralizedLangevin
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Base
)

add_library(GeneralizedLangevin
            GeneralizedLangevin.C
)

INSTALL(TARGETS GeneralizedLangevin
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
