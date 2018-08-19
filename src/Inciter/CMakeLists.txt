cmake_minimum_required(VERSION 2.8.5)

project(Inciter CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Inciter
            Transporter.C
            Partitioner.C
            FaceData.C
            Discretization.C
            Refiner.C
            Sorter.C
            MatCG.C
            DiagCG.C
            DG.C
            FluxCorrector.C
            DistFCT.C
            DiagReducer.C
            NodeDiagnostics.C
            ElemDiagnostics.C
            NodeBC.C)

target_include_directories(Inciter PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/IO
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Main
                           ${QUINOA_SOURCE_DIR}/LoadBalance
                           ${QUINOA_SOURCE_DIR}/Statistics
                           ${QUINOA_SOURCE_DIR}/LinSys
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${QUINOA_SOURCE_DIR}/PDE
                           ${QUINOA_TPL_DIR}
                           ${SEACASExodus_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${PEGTL_INCLUDE_DIRS}
                           ${HYPRE_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_DIRS}
                           ${MPI_CXX_INCLUDE_PATH}
                           ${PROJECT_BINARY_DIR}/../LinSys
                           ${PROJECT_BINARY_DIR}/../Inciter
                           ${PROJECT_BINARY_DIR}/../Main)

add_library(MeshRefinement
            AMR/mesh_adapter.C
            AMR/util.C
            AMR/edge.C
            AMR/Error.C)

target_include_directories(MeshRefinement PUBLIC
                           ${QUINOA_SOURCE_DIR}
                           ${QUINOA_SOURCE_DIR}/Base
                           ${QUINOA_SOURCE_DIR}/Control
                           ${QUINOA_SOURCE_DIR}/Mesh
                           ${QUINOA_SOURCE_DIR}/Inciter
                           ${PEGTL_INCLUDE_DIRS}
                           ${CHARM_INCLUDE_DIRS}
                           ${TPL_INCLUDE_DIR}
                           ${BRIGAND_INCLUDE_DIRS}
                           ${PROJECT_BINARY_DIR}/../Main)

addCharmModule( "transporter" "Inciter" )
addCharmModule( "partitioner" "Inciter" )
addCharmModule( "discretization" "Inciter" )
addCharmModule( "refiner" "Inciter" )
addCharmModule( "sorter" "Inciter" )
addCharmModule( "matcg" "Inciter" )
addCharmModule( "diagcg" "Inciter" )
addCharmModule( "distfct" "Inciter" )
addCharmModule( "dg" "Inciter" )

# Add extra dependency of Inciter on inciterCharmModule. This is required as
# Inciter refers to the main Charm++ proxy defined in the Charm++ module
# inciter (in Main/Inciter.C).
add_dependencies( "Inciter" "inciterCharmModule" )

## Add extra dependency of Transporter charm module on charm modules implementing
## advanced array element placement using maps. This is required since
## Transporter charm module depends on and uses advanced array element placement
## maps.
#add_dependencies( "transporterCharmModule" "linearmapCharmModule" )
#add_dependencies( "transporterCharmModule" "unsmeshmapCharmModule" )

# Add extra dependency of Inciter on Solver charm module. This is required so
# that solver.decl.h and solver.def.h are generated before Inciter including
# those.
add_dependencies( "Inciter" "solverCharmModule" )

# Add extra dependency of Transporter charm module on Partitioner charm module.
# This is required so that Partitioner.decl.h and Partitioner.def.h are
# generated before Transporter including those.
add_dependencies( "transporterCharmModule" "partitionerCharmModule" )

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
