cmake_minimum_required(VERSION 2.8.5)

project(SDE CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Particles
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/FEM
                    ${QUINOA_SOURCE_DIR}/WallTreatment
)

add_library(SDE
            Langevin.C
)

INSTALL(TARGETS SDE
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
