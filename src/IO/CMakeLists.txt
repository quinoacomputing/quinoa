cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

# Include configuration/headers/libraries common to all code
INCLUDE(${QUINOA_SOURCE_DIR}/QuinoaCommon.cmake)


include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/IO
)

add_library(IO
            IOException.C
            MeshReader.C
            MeshWriter.C
            PlotWriter.C
            GmshTxtMeshReader.C
            GmshTxtMeshWriter.C
            GmshTxtPlotWriter.C
)

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
