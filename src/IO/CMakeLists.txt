cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            Reader.C                    # readers
            GmshMeshReader.C
            STLTxtMeshReader.C
            NetgenMeshReader.C
            Writer.C                    # writers
            PDFWriter.C
            GlobWriter.C
            TxtStatWriter.C
            GmshMeshWriter.C
            SiloWriter.C
)

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
