cmake_minimum_required(VERSION 2.8.5)

project(Model CXX)

include_directories(${QUINOA_SOURCE_DIR}/Model
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Random
)

add_library(Model
            Model.C
)

INSTALL(TARGETS Model
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
