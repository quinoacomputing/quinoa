cmake_minimum_required(VERSION 2.8.5)

project(Mesh CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Mesh
)

add_library(Mesh
            Mesh.C
            StrMesh.C
            UnsMesh.C
            MeshException.C
            #Overlap.C
)

INSTALL(TARGETS Mesh
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
