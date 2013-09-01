cmake_minimum_required(VERSION 2.8.5)

project(Model CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Model/Mass
                    ${QUINOA_SOURCE_DIR}/Model/Mass/Beta
                    ${QUINOA_SOURCE_DIR}/Model/Hydro
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/SLM
                    ${QUINOA_SOURCE_DIR}/Model/Hydro/GLM
                    ${QUINOA_SOURCE_DIR}/Model/Mix
                    ${QUINOA_SOURCE_DIR}/Model/Mix/Dirichlet
                    ${QUINOA_SOURCE_DIR}/Model/Mix/GenDirichlet
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Statistics
)

add_library(Model
            Mass/Beta/Beta.C
            Hydro/SLM/SLM.C
            Hydro/GLM/GLM.C
            Mix/Dirichlet/Dirichlet.C
            Mix/GenDirichlet/GenDirichlet.C
)

INSTALL(TARGETS Model
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
