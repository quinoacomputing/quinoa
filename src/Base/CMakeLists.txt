cmake_minimum_required(VERSION 2.8.5)

project(Base CXX)

add_library(Base
            Exception.C
            Handler.C
            Timer.C
            Print.C
)

INSTALL(TARGETS Base
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
