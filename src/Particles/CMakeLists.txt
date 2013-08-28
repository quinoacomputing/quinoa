cmake_minimum_required(VERSION 2.8.5)

project(Particles CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Mesh
)

add_library(Particles
            Particles.C
)

INSTALL(TARGETS Particles
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
