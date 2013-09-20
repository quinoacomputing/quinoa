cmake_minimum_required(VERSION 2.8.5)

project(Model CXX)

add_library(Model
            Mass/Beta/Beta.C
            Hydro/SLM/SLM.C
            Hydro/GLM/GLM.C
            Mix/Dirichlet/Dirichlet.C
            Mix/GenDirichlet/GenDirichlet.C
)

INSTALL(TARGETS Model
  	RUNTIME DESTINATION bin COMPONENT Runtime
  	LIBRARY DESTINATION lib COMPONENT Runtime
 	ARCHIVE DESTINATION lib COMPONENT Development
)
