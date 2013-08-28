cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Mesh
)

add_library(Mesh
            GmshMesh.C
            STLMesh.C
)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
