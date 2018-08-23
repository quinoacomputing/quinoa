cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            PDFWriter.C
            TxtStatWriter.C
            DiagWriter.C)

target_include_directories(IO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${QUINOA_DIRS}
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(IO PROPERTIES LIBRARY_OUTPUT_NAME quinoa_io)

INSTALL(TARGETS IO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

# Native IO libraries do not require third-party libraries
add_library(NativeMeshIO
            MeshFactory.C
            GmshMeshReader.C
            STLTxtMeshReader.C
            NetgenMeshReader.C
            ASCMeshReader.C
            GmshMeshWriter.C
            NetgenMeshWriter.C
            #SiloWriter.C
)

target_include_directories(NativeMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Control
                           ${TPL_INCLUDE_DIR}
                           ${CHARM_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PROJECT_BINARY_DIR}/../Main)

set_target_properties(NativeMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_nativemeshio)

INSTALL(TARGETS NativeMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(ExodusIIMeshIO
            ExodusIIMeshReader.C
            ExodusIIMeshWriter.C)

target_include_directories(ExodusIIMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH})

set_target_properties(ExodusIIMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_exodusiimeshio)

INSTALL(TARGETS ExodusIIMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

if (HAS_OMEGA_H)
  add_library(OmegaHMeshIO
              Omega_h_MeshReader.C)

  target_include_directories(OmegaHMeshIO PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Mesh
                             ${OMEGA_H_INCLUDE_DIRS}
                             ${MPI_CXX_INCLUDE_DIRS}
                             ${MPI_CXX_INCLUDE_PATH})

  set_target_properties(OmegaHMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                        quinoa_omegahmeshio)
  
  INSTALL(TARGETS OmegaHMeshIO
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
  )
endif()

if (HAS_ROOT)
  add_library(RootMeshIO
              FileConvWriter.C
              RootMeshWriter.C)
  
  target_include_directories(RootMeshIO PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Mesh
                             ${SEACASExodus_INCLUDE_DIRS}
                             ${ROOT_INCLUDE_DIRS}
                             ${MPI_CXX_INCLUDE_DIRS}
                             ${MPI_CXX_INCLUDE_PATH})

  set_target_properties(RootMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                        quinoa_rootmeshio)
  
  INSTALL(TARGETS RootMeshIO
          RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
          LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
          ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
  )
endif()

add_library(HyperMeshIO
            HyperMeshReader.C)

target_include_directories(HyperMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${TPL_INCLUDE_DIR}
                           ${PUGIXML_INCLUDE_DIRS})

set_target_properties(HyperMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_hypermeshio)

INSTALL(TARGETS HyperMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(MeshDetect
            MeshDetect.C)

target_include_directories(MeshDetect PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base)

set_target_properties(MeshDetect PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_meshdetect)

INSTALL(TARGETS MeshDetect
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(ParticleIO
            H5PartWriter.C)

target_include_directories(ParticleIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${H5PART_INCLUDE_DIRS}
                           ${HDF5_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH})

set_target_properties(ParticleIO PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_particleio)

INSTALL(TARGETS ParticleIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
