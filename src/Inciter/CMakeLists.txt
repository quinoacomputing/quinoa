project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Transporter.cpp
            Partitioner.cpp
            FaceData.cpp
            Discretization.cpp
            Refiner.cpp
            Sorter.cpp
            DiagCG.cpp
            ALECG.cpp
            ALE.cpp
            OversetFE.cpp
            Ghosts.cpp
            DG.cpp
            FV.cpp
            FluxCorrector.cpp
            DistFCT.cpp
            DiagReducer.cpp
            HoleReducer.cpp
            NodeDiagnostics.cpp
            ElemDiagnostics.cpp
            FieldOutput.cpp
            NodeBC.cpp)

target_include_directories(Inciter PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/IO
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Main
                           ${QUINOA_SOURCE_DIR}/LoadBalance
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${QUINOA_SOURCE_DIR}/LinearSolver
                           ${QUINOA_SOURCE_DIR}/Transfer
                           ${QUINOA_TPL_DIR}
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${NETCDF_INCLUDES}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HYPRE_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${EXAM2M_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${CBLAS_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../LinearSolver
                           ${PROJECT_BINARY_DIR}/../Base
                           ${PROJECT_BINARY_DIR}/../IO
                           ${PROJECT_BINARY_DIR}/../Mesh
                           ${PROJECT_BINARY_DIR}/../Transfer
                           ${PROJECT_BINARY_DIR}/../Main)

add_library(MeshRefinement
            PUPAMR.cpp
            AMR/mesh_adapter.cpp
            AMR/util.cpp
            AMR/edge.cpp
            AMR/Error.cpp)

target_include_directories(MeshRefinement PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${HIGHWAYHASH_INCLUDE_DIRS}
                           ${TPL_INCLUDE_DIR}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${CBLAS_INCLUDE_DIRS}
                           ${LAPACKE_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

addCharmModule( "transporter" "Inciter" )
addCharmModule( "partitioner" "Inciter" )
addCharmModule( "discretization" "Inciter" )
addCharmModule( "refiner" "Inciter" )
addCharmModule( "sorter" "Inciter" )
addCharmModule( "diagcg" "Inciter" )
addCharmModule( "alecg" "Inciter" )
addCharmModule( "oversetfe" "Inciter" )
addCharmModule( "distfct" "Inciter" )
addCharmModule( "dg" "Inciter" )
addCharmModule( "fv" "Inciter" )
addCharmModule( "ale" "Inciter" )
addCharmModule( "ghosts" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.cpp).
add_dependencies( "Inciter" "inciterCharmModule" )

## Add extra dependency of Transporter charm module on charm modules implementing
## advanced array element placement using maps. This is required since
## Transporter charm module depends on and uses advanced array element placement
## maps.
#add_dependencies( "transporterCharmModule" "linearmapCharmModule" )
#add_dependencies( "transporterCharmModule" "unsmeshmapCharmModule" )

# Add extra dependency of Transporter charm module on Partitioner charm module.
# This is required so that partitioner.decl.h and partitioner.def.h are
# generated before Transporter including those.
add_dependencies( "transporterCharmModule" "partitionerCharmModule" )
add_dependencies( "transporterCharmModule" "m2mtransferCharmModule" )

# Add extra dependency of Discretization charm module on MeshWriter charm
# module. This is so meshwriter.decl.h and meshwriter.def.h are generated
# before Discretization which needs meshwriter's type information.
add_dependencies( "discretizationCharmModule" "meshwriterCharmModule" )

set_target_properties(Inciter PROPERTIES LIBRARY_OUTPUT_NAME quinoa_inciter)
set_target_properties(MeshRefinement PROPERTIES LIBRARY_OUTPUT_NAME
                      quinoa_mesh_refinement)

INSTALL(TARGETS Inciter
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

INSTALL(TARGETS MeshRefinement
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
