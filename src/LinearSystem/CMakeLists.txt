cmake_minimum_required(VERSION 2.8.5)

project(LinearSystem CXX)

add_library(LinearSystem
            ZoltanInterOp.C
            HypreInterOp.C
)

removeWarnings("ZoltanInterOp.C;HypreInterOp.C")

INSTALL(TARGETS LinearSystem
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
