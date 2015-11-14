

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(phdmesh_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/phdmesh_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include)

#
# Element files
#

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include/element)

SET(ELEMENT_HEADERS ${HEADERS}
  ../include/element/Basic_Topologies.hpp
  ../include/element/CellTopology.hpp
  ../include/element/Declarations.hpp
  ../include/element/Dimensions.hpp
  ../include/element/Hexahedron_Topologies.hpp
  ../include/element/Pyramid_Topologies.hpp
  ../include/element/Quadrilateral_Topologies.hpp
  ../include/element/Stencils.hpp
  ../include/element/Tetrahedron_Topologies.hpp
  ../include/element/Triangle_Topologies.hpp
  ../include/element/Wedge_Topologies.hpp
  )

SET(ELEMENT_SOURCES
  element/ElementDeclarations.cpp
  element/ElementDimensions.cpp
  element/ElementTopologies.cpp
  )

#
# Mesh files
#

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../include/mesh)

SET(MESH_HEADERS ${HEADERS}
  ../include/mesh/Assemble.hpp
  ../include/mesh/BulkData.hpp
  ../include/mesh/Comm.hpp
  ../include/mesh/Entity.hpp
  ../include/mesh/EntityComm.hpp
  ../include/mesh/Field.hpp
  ../include/mesh/FieldDeclarations.hpp
  ../include/mesh/FieldData.hpp
  ../include/mesh/FieldParallel.hpp
  ../include/mesh/FieldTraits.hpp
  ../include/mesh/Gather.hpp
  ../include/mesh/Kernel.hpp
  ../include/mesh/MetaData.hpp
  ../include/mesh/Part.hpp
  ../include/mesh/Proximity.hpp
  ../include/mesh/Relation.hpp
  ../include/mesh/Types.hpp
  )

SET(MESH_SOURCES
  mesh/MeshAssemble.cpp
  mesh/MeshAura.cpp
  mesh/MeshBulkData.cpp
  mesh/MeshComm.cpp
  mesh/MeshEntity.cpp
  mesh/MeshEntityComm.cpp
  mesh/MeshField.cpp
  mesh/MeshFieldParallel.cpp
  mesh/MeshFieldTraits.cpp
  mesh/MeshKernel.cpp
  mesh/MeshMetaData.cpp
  mesh/MeshPart.cpp
  mesh/MeshProximity.cpp
  mesh/MeshRebalance.cpp
  mesh/MeshRelation.cpp
  mesh/MeshSharing.cpp
  mesh/MeshTypes.cpp
  )

#
# Util files
#

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/util)

SET(UTIL_HEADERS ${HEADERS}
  ../include/util/Array.hpp
  ../include/util/ArrayPrivate.hpp
  ../include/util/Basics.hpp
  ../include/util/CSet.hpp
  ../include/util/FixedPoolAlloc.hpp
  ../include/util/IndexList.hpp
  ../include/util/NamedValue.hpp
  ../include/util/NumericEnum.hpp
  ../include/util/OctTree.hpp
  ../include/util/OctTreeOps.hpp
  ../include/util/PairIter.hpp
  ../include/util/Parallel.hpp
  ../include/util/ParallelComm.hpp
  ../include/util/ParallelIndex.hpp
  ../include/util/ParallelInputStream.hpp
  ../include/util/ParallelReduce.hpp
  ../include/util/Setv.hpp
  ../include/util/TestDriver.hpp
  ../include/util/TPI.h
  ../include/util/TPI.hpp
  ../include/util/TypeList.hpp
  ../include/util/TypeName.hpp
  )

SET(UTIL_SOURCES ${HEADERS}
  util/ArrayPrivate.cpp
  util/CSet.cpp
  util/FixedPoolAlloc.cpp
  util/NamedValue.cpp
  util/NumericEnum.cpp
  util/OctTree.cpp
  util/OctTreeOps.cpp
  util/ParallelComm.cpp
  util/Parallel.cpp
  util/ParallelIndex.cpp
  util/ParallelInputStream.cpp
  util/ParallelReduce.cpp
  util/Setv.cpp
  util/TestDriver.cpp
  util/TPI_Concurrency.c
  util/TPI.c
  util/TPI_Walltime.c
  util/TypeName.cpp
  )

#
# yaml files
#

SET(YAML_HEADERS ${HEADERS}
  yaml/YAML_Doc.hpp
  yaml/YAML_Element.hpp
  )

SET(YAML_SOURCES ${HEADERS}
  yaml/YAML_Doc.cpp
  yaml/YAML_Element.cpp
  )

#does this need to be enclosed by the exodus option?

#
# Mesh IO files
#

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/mesh_io)

SET(MESH_IO_HEADERS ${HEADERS}
  ../include/mesh_io/ExoII.hpp
  )

SET(MESH_IO_SOURCES ${HEADERS}
  mesh_io/MeshIO_ExoII.cpp
  )

#
# TX Blas files
#

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/txblas)

SET(TXBLAS_HEADERS ${HEADERS}
  ../include/txblas/CR4Matrix.hpp
  ../include/txblas/cr4_mxv.h
  ../include/txblas/reduction.h
  ../include/txblas/Reduction.hpp
  )

SET(TXBLAS_SOURCES ${HEADERS}
  txblas/txblas_axpby.c
  txblas/txblas_CR4Matrix.cpp
  txblas/txblas_cr4_mxv.c
  txblas/txblas_reduction.c
  )

IF (TPL_ENABLE_MPI)
ENDIF()


#
# C) Define the targets for package's library(s)
#
TRIBITS_ADD_LIBRARY(
  phdmesh_util
  HEADERS ${UTIL_HEADERS}
  SOURCES ${UTIL_SOURCES}
  )

TRIBITS_ADD_LIBRARY(
  phdmesh_mesh
  HEADERS ${MESH_HEADERS}
  SOURCES ${MESH_SOURCES}
  DEPLIBS phdmesh_util
  )

TRIBITS_ADD_LIBRARY(
  phdelement_mesh
  HEADERS ${ELEMENT_HEADERS}
  SOURCES ${ELEMENT_SOURCES}
  DEPLIBS phdmesh_mesh
  )

TRIBITS_ADD_LIBRARY(
  phdmesh_yaml
  HEADERS ${YAML_HEADERS}
  SOURCES ${YAML_SOURCES}
  )

TRIBITS_ADD_LIBRARY(
  phdmesh_meshio
  HEADERS ${MESH_IO_HEADERS}
  SOURCES ${MESH_IO_SOURCES}
  DEPLIBS phdmesh_mesh phdelement_mesh
  )

TRIBITS_ADD_LIBRARY(
  phdmesh_txblas
  HEADERS ${TXBLAS_HEADERS}
  SOURCES ${TXBLAS_SOURCES}
  DEPLIBS phdmesh_util
  )
