cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            PDFWriter.C
            TxtStatWriter.C
            DiagWriter.C
            H5PartWriter.C
)

INSTALL(TARGETS IO
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

add_library(MeshIO
            MeshFactory.C
            GmshMeshReader.C
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ExodusIIMeshReader.C
            HyperMeshReader.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            ExodusIIMeshWriter.C
            #SiloWriter.C
)

INSTALL(TARGETS MeshIO
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)

# Include function for adding Charm++ modules
include(charm)

add_library(AsyncIO
            ParticleWriter.C
)

addCharmModule( "particlewriter" "AsyncIO" )

# Add extra dependency of AsyncIO on Transporter charm module. This is required
# so that transporter.decl.h and transporter.def.h are generated before AsyncIO
# including those.
add_dependencies( "AsyncIO" "transporterCharmModule" )

INSTALL(TARGETS AsyncIO
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib COMPONENT Development
)
