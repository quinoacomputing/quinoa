cmake_minimum_required(VERSION 2.8.5)

project(LinSys CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(LinSys
            DiagReducer.C
            Solver.C)

addCharmModule( "solver" "LinSys" )

# Add extra dependency of LinSys on charm modules implementing advanced array
# element placement using maps. This is required to pull in the dependencies of
# the possible specializations of LinSys using various map placement strategies.
# For example, LinSys is specialized with Transporter which in turn depends on
# array placement maps.
add_dependencies( "LinSys" "linearmapCharmModule" )
add_dependencies( "LinSys" "unsmeshmapCharmModule" )

# Add extra dependency of Solver charm module on Transporter charm module.  This
# is required to make sure transporter.decl.h is generated before solver.ci is
# processed, as Solver is specialized on CProxy_Transporter.
add_dependencies( "solverCharmModule" "transporterCharmModule" )

# Add extra dependency of Solver charm module on CG charm module.  This is
# required to make sure matcg.decl.h is generated before solver.ci is processed,
# as Solver is specialized on CProxy_MatCG.
add_dependencies( "solverCharmModule" "matcgCharmModule" )

set_target_properties(LinSys PROPERTIES LIBRARY_OUTPUT_NAME quinoa_linsys)

INSTALL(TARGETS LinSys
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
