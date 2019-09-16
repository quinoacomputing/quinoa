cmake_minimum_required(VERSION 2.8.5)

project(IO CXX)

add_library(IO
            PDFWriter.cpp
            TxtStatWriter.cpp
            DiagWriter.cpp)

target_include_directories(IO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${NETCDF_INCLUDES}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
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
            MeshFactory.cpp
            GmshMeshReader.cpp
            STLTxtMeshReader.cpp
            NetgenMeshReader.cpp
            ASCMeshReader.cpp
            GmshMeshWriter.cpp
            NetgenMeshWriter.cpp
            #SiloWriter.cpp
)

target_include_directories(NativeMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Control
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${NETCDF_INCLUDES}
                           ${TPL_INCLUDE_DIR}
                           ${CHARM_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
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
            ExodusIIMeshReader.cpp
            ExodusIIMeshWriter.cpp)

target_include_directories(ExodusIIMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${NETCDF_INCLUDES}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
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
              Omega_h_MeshReader.cpp)

  target_include_directories(OmegaHMeshIO PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Mesh
                             ${TPL_INCLUDE_DIR}
                             ${OMEGA_H_INCLUDE_DIRS}
                             ${HIGHWAYHASH_INCLUDE_DIRS}
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
              FileConvWriter.cpp
              RootMeshWriter.cpp)
  
  target_include_directories(RootMeshIO PUBLIC
                             ${QUINOA_SOURCE_DIR}
                             ${QUINOA_SOURCE_DIR}/Base
                             ${QUINOA_SOURCE_DIR}/Mesh
                             ${SEACASExodus_INCLUDE_DIRS}
                             ${NETCDF_INCLUDES}
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
            HyperMeshReader.cpp)

target_include_directories(HyperMeshIO PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${TPL_INCLUDE_DIR}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${PUGIXML_INCLUDE_DIRS})

set_target_properties(HyperMeshIO PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_hypermeshio)

INSTALL(TARGETS HyperMeshIO
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

add_library(MeshDetect
            MeshDetect.cpp)

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

add_library(MeshWriter
            MeshWriter.cpp)

target_include_directories(MeshWriter PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Control
                           ${PROJECT_BINARY_DIR}/../Main
                           ${PROJECT_BINARY_DIR}/../IO
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${NETCDF_INCLUDES}
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${ROOT_INCLUDE_DIRS}
                           ${BRIGAND_INCLUDE_DIRS})

set_target_properties(MeshWriter PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_meshwriter)

INSTALL(TARGETS MeshWriter
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

addCharmModule( "meshwriter" "MeshWriter" )

# The following lines are commented out so that meshconv can be cuilt without
# h5part. Once ParticleIO is hooked in again, these lines will be uncommented.
#add_library(ParticleIO
#            H5PartWriter.cpp)
#
#target_include_directories(ParticleIO PUBLIC
#                           ${QUINOA_SOURCE_DIR}
#                           ${QUINOA_SOURCE_DIR}/Base
#                           ${H5PART_INCLUDE_DIRS}
#                           ${HDF5_INCLUDE_DIRS}
#                           ${MPI_CXX_INCLUDE_DIRS}
#                           ${MPI_CXX_INCLUDE_PATH})
#
#set_target_properties(ParticleIO PROPERTIES LIBRARY_OUTPUT_NAME
#                      quinoa_particleio)
#
#INSTALL(TARGETS ParticleIO
#        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
#        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
#        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
#)
