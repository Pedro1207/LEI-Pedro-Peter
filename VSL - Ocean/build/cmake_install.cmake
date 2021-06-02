# Install script for directory: F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/bin")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/MinGW/bin/objdump.exe")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/contrib/glew/cmake_install.cmake")
  include("F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/contrib/tinyxml/cmake_install.cmake")
  include("F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/contrib/freeglut-3.0.0/cmake_install.cmake")
  include("F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/VSL/cmake_install.cmake")
  include("F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/demo/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "F:/Main/Uni/LEI - GIT/LEI-Pedro-Peter/VSL - Ocean/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
