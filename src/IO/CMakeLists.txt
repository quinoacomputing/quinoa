cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            PDFWriter.C
            TxtStatWriter.C
            DiagWriter.C
            H5PartWriter.C
)

set_target_properties(IO PROPERTIES LIBRARY_OUTPUT_NAME quinoa_io)

INSTALL(TARGETS IO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

if (HAS_ROOT)
  set(ROOT_WRITER RootMeshWriter.C)
  set(FILE_CONVERTER FileConvWriter.C)
endif()

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
	    ${ROOT_WRITER}
	    ${FILE_CONVERTER}
            #SiloWriter.C
)

set_target_properties(MeshIO PROPERTIES LIBRARY_OUTPUT_NAME quinoa_meshio)

INSTALL(TARGETS MeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
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

set_target_properties(AsyncIO PROPERTIES LIBRARY_OUTPUT_NAME quinoa_asyncio)

INSTALL(TARGETS AsyncIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
