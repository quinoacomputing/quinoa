cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

include_directories(${QUINOA_SOURCE_DIR}/Base
                    ${QUINOA_SOURCE_DIR}/Mesh
                    ${QUINOA_SOURCE_DIR}/IO
                    ${QUINOA_SOURCE_DIR}/Control
                    ${QUINOA_SOURCE_DIR}/Statistics
                    ${TPL_DIR}/include
)

add_library(IO
            Reader.C                    # readers
            GmshTxtMeshReader.C
            STLTxtMeshReader.C
            Writer.C                    # writers
            PDFWriter.C
            GlobWriter.C
            TxtStatWriter.C
            GmshTxtPlotWriter.C
            GmshTxtMeshWriter.C
            SiloWriter.C
)

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
