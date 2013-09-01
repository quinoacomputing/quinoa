cmake_minimum_required(VERSION 2.8.5)

project(Physics CXX)

include_directories(${QUINOA_SOURCE_DIR}/Physics
                    ${QUINOA_SOURCE_DIR}/Physics/HomMix
                    ${QUINOA_SOURCE_DIR}/Physics/HomHydro
                    ${QUINOA_SOURCE_DIR}/Physics/HomRT
                    ${QUINOA_SOURCE_DIR}/Physics/SPINSFlow
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Paradigm
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model/Mass/Beta
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model/Mix/Dirichlet
                    ${QUINOA_SOURCE_DIR}/Model/Mix/GenDirichlet
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/SLM
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/GLM
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Random
)

add_library(Physics
            Physics.C
            HomMix/HomMix.C
            HomHydro/HomHydro.C
            HomRT/HomRT.C
            SPINSFlow/SPINSFlow.C
)

INSTALL(TARGETS Physics
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
