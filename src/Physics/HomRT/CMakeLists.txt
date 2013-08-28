cmake_minimum_required(VERSION 2.8.5)

project(HomRT CXX)

include_directories(${QUINOA_SOURCE_DIR}/Physics/HomRT
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model/Mass/Beta
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/SimplifiedLangevin
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/GeneralizedLangevin
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Control
)

add_library(HomRT
            HomRT.C
)

INSTALL(TARGETS HomRT
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
