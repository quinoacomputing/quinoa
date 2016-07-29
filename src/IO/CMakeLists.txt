cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(IO
            MeshFactory.C
            GmshMeshReader.C            # readers
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ExodusIIMeshReader.C
            HyperMeshReader.C
            PDFWriter.C                 # writers
            TxtStatWriter.C
            DiagWriter.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            ExodusIIMeshWriter.C
            H5PartWriter.C
            ParticleWriter.C
#            SiloWriter.C
)

addCharmModule( "particlewriter" "IO" )

INSTALL(TARGETS IO
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
