cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            GmshMeshReader.C            # readers
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ExodusIIMeshReader.C
            PDFWriter.C                 # writers
            GlobWriter.C
            TxtStatWriter.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            ExodusIIMeshWriter.C
            SiloWriter.C
)

removeWarnings( "TxtStatWriter.C" )

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
