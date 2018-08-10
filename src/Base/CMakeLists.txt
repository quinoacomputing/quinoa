cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

# Include function for adding Charm++ modules
include(charm)

add_library(Base
            Exception.C
            ProcessException.C
            ProcessControl.C
            LoadDistributor.C
            Timer.C
            Reader.C
            Writer.C
            VectorReducer.C
            Table.C
            ChareStateCollector.C
)

addCharmModule( "charestatecollector" "Base" )

set_target_properties(Base PROPERTIES LIBRARY_OUTPUT_NAME quinoa_base)

INSTALL(TARGETS Base
       RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
       LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
       ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
