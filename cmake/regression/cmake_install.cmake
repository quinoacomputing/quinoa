# Install script for directory: /home/aditya/quinoa/regression

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/aditya/quinoa/cmake/regression/walker/Dirichlet/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/GeneralizedDirichlet/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/SkewNormal/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/DiagOrnsteinUhlenbeck/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/OrnsteinUhlenbeck/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/Beta/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/Gamma/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/NumFracBeta/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/MassFracBeta/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/walker/MixMassFracBeta/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/rngtest/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/meshconv/gmsh_output/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/meshconv/netgen_output/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/meshconv/exo_output/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/inciter/asynclogic/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/inciter/fct/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/inciter/cfl/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/inciter/transport/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/regression/inciter/amr/initial/uniform/cmake_install.cmake")

endif()

