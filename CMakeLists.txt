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

list(
  APPEND _py_dims 1 2 3 D
    )

list(
  APPEND _py_metrics euclidean homogeneous conformal minkowski
    )

foreach(_metric IN ITEMS ${_py_metrics})
  foreach(_dim IN ITEMS ${_py_dims})
    message(Building\ ${_metric}${_dim})
    python_add_module(${_metric}${_dim} ./python/src/${_metric}/${_metric}${_dim}.cpp)
    target_link_libraries(${_metric}${_dim} PUBLIC Boost::${_boost_component_found} ${PYTHON_LIBRARIES})
    target_include_directories(${_metric}${_dim} PRIVATE ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
    set_target_properties(${_metric}${_dim} 
      PROPERTIES
      PREFIX ""
      ARCHIVE_OUTPUT_DIRECTORY "tbgal/${_metric}"
      LIBRARY_OUTPUT_DIRECTORY "tbgal/${_metric}"
      RUNTIME_OUTPUT_DIRECTORY "tbgal/${_metric}"
    )
  add_custom_command(TARGET ${_metric}${_dim} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory tbgal)
  add_custom_command(TARGET ${_metric}${_dim} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch tbgal/__init__.py)
  endforeach()
endforeach()


message(Building\ SignedPQ)
python_add_module(PQ ./python/src/signed/PQ.cpp)
target_link_libraries(PQ PUBLIC Boost::${_boost_component_found} ${PYTHON_LIBRARIES})
target_include_directories(PQ PRIVATE ${PYTHON_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
set_target_properties(PQ 
  PROPERTIES
  PREFIX ""
  ARCHIVE_OUTPUT_DIRECTORY "tbgal/signed"
  LIBRARY_OUTPUT_DIRECTORY "tbgal/signed"
  RUNTIME_OUTPUT_DIRECTORY "tbgal/signed"
)
add_custom_command(TARGET PQ PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory tbgal)
add_custom_command(TARGET PQ PRE_BUILD COMMAND ${CMAKE_COMMAND} -E touch tbgal/__init__.py)


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
