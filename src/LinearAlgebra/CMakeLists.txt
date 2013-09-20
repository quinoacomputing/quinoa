cmake_minimum_required(VERSION 2.8.5)

project(LinearAlgebra CXX)

add_library(LinearAlgebra
            #CG.C
            #Matrix3.C
            SymCompRowMatrix.C
)

INSTALL(TARGETS LinearAlgebra
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
