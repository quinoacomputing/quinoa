cmake_minimum_required(VERSION 2.8.5)

project(Particles CXX)

add_library(Particles
            Tracker.C
)

INSTALL(TARGETS Particles
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)

# Add extra dependency of Particles on ParticleWriter charm module. This is
# required so that particlewriter.decl.h and particlewriter.def.h are generated
# before Tracker including those.
add_dependencies( "Particles" "particlewriterCharmModule" )

# Add extra dependency of Particles on Transporter charm module. This is
# required so that transporter.decl.h and transporter.def.h are generated before
# Tracker including those.
add_dependencies( "Particles" "transporterCharmModule" )
