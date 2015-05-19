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
            LinearMap.C
#            Msg.C
)

#addCharmModule( "msg" "Base" )
addCharmModule( "linearmap" "Base" )
removeWarnings( "ProcessException.C;Timer.C;LinearMap.C")

INSTALL(TARGETS Base
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
