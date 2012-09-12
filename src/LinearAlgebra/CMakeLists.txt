cmake_minimum_required(VERSION 2.8.5)

project(LinearAlgebra CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/LinearAlgebra
)

add_library(LinearAlgebra
            #CG.C
            #Matrix3.C
            SparseMatrix.C  
            SymCompressedRowMatrix.C
)

INSTALL(TARGETS LinearAlgebra
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
