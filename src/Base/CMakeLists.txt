cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

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
#            Msg.C
)

#addCharmModule( "msg" "Base" )

INSTALL(TARGETS Base
  	RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT Runtime
  	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Runtime
 	ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Development
)
