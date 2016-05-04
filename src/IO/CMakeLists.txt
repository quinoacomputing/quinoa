cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            MeshFactory.C
            GmshMeshReader.C            # readers
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ExodusIIMeshReader.C
            HyperMeshReader.C
            PDFWriter.C                 # writers
            TxtStatWriter.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            ExodusIIMeshWriter.C
#            SiloWriter.C
)

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
