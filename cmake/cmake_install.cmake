# Install script for directory: /home/aditya/quinoa/src

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
  include("/home/aditya/quinoa/cmake/Main/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/UnitTest/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Base/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Control/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/DiffEq/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/PDE/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Walker/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Inciter/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/LoadBalance/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/LinSys/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/IO/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Mesh/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/RNG/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/RNGTest/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Statistics/cmake_install.cmake")
  include("/home/aditya/quinoa/cmake/Particles/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/aditya/quinoa/cmake/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
