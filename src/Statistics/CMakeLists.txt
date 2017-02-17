cmake_minimum_required(VERSION 2.8.5)

project(Statistics CXX)

add_library(Statistics
            Statistics.C
            PDFReducer.C
)

set_target_properties(Statistics PROPERTIES
                      LIBRARY_OUTPUT_NAME quinoa_statistics)

INSTALL(TARGETS Statistics
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
