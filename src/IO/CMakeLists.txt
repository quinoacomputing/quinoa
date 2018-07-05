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

if (HAS_OMEGA_H)
  set(OMEGA_H_READER Omega_h_MeshReader.C)
endif()

add_library(MeshIO
            MeshFactory.C
            GmshMeshReader.C
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ExodusIIMeshReader.C
            HyperMeshReader.C
            ASCMeshReader.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            ExodusIIMeshWriter.C
	    ${ROOT_WRITER}
	    ${FILE_CONVERTER}
            ${OMEGA_H_READER}
            #SiloWriter.C
)

set_target_properties(MeshIO PROPERTIES LIBRARY_OUTPUT_NAME quinoa_meshio)

INSTALL(TARGETS MeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
