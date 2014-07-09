cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

include(charm)

add_library(Base
            Exception.C
            Handler.C
            Timer.C
            RNGPrint.C
            QuinoaPrint.C
#            Msg.C
)

add_dependencies(Base testu01suiteCharmModule)
add_dependencies(Base testu01CharmModule)

#addCharmModule( "msg" "Base" )
removeWarnings( "QuinoaPrint.C;RNGPrint.C;Handler.C" )

INSTALL(TARGETS Base
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
