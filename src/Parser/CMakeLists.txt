cmake_minimum_required(VERSION 2.8.5)

project(Parser CXX)

include_directories(${QUINOA_SOURCE_DIR}/Parser
                    ${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Geometry
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/IO
                    ${TPL_DIR}/include
)

add_library(Parser
            Parser.C
)

INSTALL(TARGETS Parser
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
