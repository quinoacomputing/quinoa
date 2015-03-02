cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

include(charm)

add_library(Base
            Exception.C
            ProcessException.C
            LoadDistributor.C
            Timer.C
            Reader.C
            Writer.C
#            Msg.C
)

#addCharmModule( "msg" "Base" )
removeWarnings( "LoadDistributor.C;ProcessException.C")

INSTALL(TARGETS Base
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
