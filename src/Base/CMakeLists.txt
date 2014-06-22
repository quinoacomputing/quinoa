cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

include(charm)

add_library(Base
            Exception.C
            Handler.C
            Timer.C
            Print.C
            RNGPrint.C
            QuinoaPrint.C
            Msg.C
)

addCharmModule( "msg" "Base" )
removeWarnings( "QuinoaPrint.C;RNGPrint.C;Handler.C;Msg.C" )

INSTALL(TARGETS Base
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
