# Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
# 
# authors    : Sousa, Eduardo V.
#              Fernandes, Leandro A. F.
# repository : https://github.com/Prograf-UFF/TbGAL
# 
# This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
# 
# TbGAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# TbGAL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with TbGAL. If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.10)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 20200123)

project(TbGAL
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
  DESCRIPTION "TbGAL: Tensor-based Geometric Algebra Library"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 NO_MODULE REQUIRED)

# for testing we will need the python interpreter
find_package(PythonInterp REQUIRED)

# we require python development headers
find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EXACT REQUIRED)

#file(MAKE_DIRECTORY reports)
#file(MAKE_DIRECTORY tbgal)
#file(WRITE tbgal/__init.py__)

# now search for the boost component
# depending on the boost version it is called either python,
# python2, python27, python3, python36, python37, ...

list(
  APPEND _components
    python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
    python${PYTHON_VERSION_MAJOR}
    python
  )

set(_boost_component_found "")

foreach(_component IN ITEMS ${_components})
  find_package(Boost COMPONENTS ${_component})
  if(Boost_FOUND)
    set(_boost_component_found ${_component})
    break()
  endif()
endforeach()

if(_boost_component_found STREQUAL "")
  message(FATAL_ERROR "No matching Boost.Python component found")
endif()

python_add_module(conformal ./python/src/conformal.cpp)
python_add_module(euclidean ./python/src/euclidean.cpp)
python_add_module(homogeneous ./python/src/homogeneous.cpp)

target_link_libraries(conformal PUBLIC Boost::${_boost_component_found} ${PYTHON_LIBRARIES})
target_link_libraries(euclidean PUBLIC Boost::${_boost_component_found} ${PYTHON_LIBRARIES})
target_link_libraries(homogeneous PUBLIC Boost::${_boost_component_found} ${PYTHON_LIBRARIES})

target_include_directories(conformal PRIVATE ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
target_include_directories(euclidean PRIVATE ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
target_include_directories(homogeneous PRIVATE ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})

set_target_properties(conformal
  PROPERTIES
    PREFIX ""
    ARCHIVE_OUTPUT_DIRECTORY "tbgal"
    LIBRARY_OUTPUT_DIRECTORY "tbgal"
    RUNTIME_OUTPUT_DIRECTORY "tbgal"
  )

set_target_properties(homogeneous
  PROPERTIES
    PREFIX ""
    ARCHIVE_OUTPUT_DIRECTORY "tbgal"
    LIBRARY_OUTPUT_DIRECTORY "tbgal"
    RUNTIME_OUTPUT_DIRECTORY "tbgal"
  )

set_target_properties(euclidean
  PROPERTIES
    PREFIX ""
    ARCHIVE_OUTPUT_DIRECTORY "tbgal"
    LIBRARY_OUTPUT_DIRECTORY "tbgal"
    RUNTIME_OUTPUT_DIRECTORY "tbgal"
  )


add_custom_command(TARGET conformal PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory tbgal)
add_custom_command(TARGET conformal PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch tbgal/__init.py__)
add_custom_command(TARGET homogeneous PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory tbgal)
add_custom_command(TARGET homogeneous PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch tbgal/__init.py__)
add_custom_command(TARGET euclidean PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory tbgal)
add_custom_command(TARGET euclidean PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch tbgal/__init.py__)

if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(conformal PRIVATE "-Wall -Wextra")
    target_compile_options(homogeneous PRIVATE "-Wall -Wextra")
    target_compile_options(euclidean PRIVATE "-Wall -Wextra")
endif()
if ( MSVC )
    target_compile_options(conformal PRIVATE "/W4")
    target_compile_options(homogeneous PRIVATE "/W4")
    target_compile_options(euclidean PRIVATE "/W4")
endif()

message(STATUS "Eigen3_DIR = ${EIGEN3_INCLUDE_DIR}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")

file(INSTALL "./cpp/include" DESTINATION "${CMAKE_INSTALL_PREFIX}" FILES_MATCHING PATTERN "*.hpp")

configure_file("./cmake/TbGALConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")

configure_file("./cmake/TbGALConfigVersion.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")
