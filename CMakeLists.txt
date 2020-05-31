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

cmake_minimum_required(VERSION 3.14)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 20200531)

project(TbGAL
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
  DESCRIPTION "TbGAL: Tensor-Based Geometric Algebra Library"
  HOMEPAGE_URL https://github.com/Prograf-UFF/TbGAL
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 NO_MODULE REQUIRED)

message(STATUS "Installing C++ front-end")

file(INSTALL ./cpp/include DESTINATION ${CMAKE_INSTALL_PREFIX} FILES_MATCHING PATTERN "*.hpp")

configure_file(./cmake/TbGALConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake @ONLY NEWLINE_STYLE UNIX)
file(INSTALL ${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/tbgal)

configure_file(./cmake/TbGALConfigVersion.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake @ONLY NEWLINE_STYLE UNIX)
file(INSTALL ${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/tbgal)

message(STATUS "Installing C++ front-end - done")

message(STATUS "Configuring Python front-end")

find_package(Python COMPONENTS Interpreter Development NumPy)
if(NOT Python_FOUND)
  message(STATUS "Configuring Python front-end - error (Python.Interpreter, Python.Developmemnt, or Python.NumPy components found)")
else()
  if(WIN32)
    set(Boost_USE_MULTITHREADED ON)
    set(Boost_USE_STATIC_LIBS OFF)
    add_definitions(-DBOOST_ALL_NO_LIB)
    add_definitions(-DBOOST_ALL_DYN_LINK)
  endif()

  foreach(_boost_python_suffix IN ITEMS "${Python_VERSION_MAJOR}${Python_VERSION_MINOR}" "${Python_VERSION_MAJOR}" "")
    find_package(Boost COMPONENTS python${_boost_python_suffix} numpy${_boost_python_suffix})
    if(Boost_FOUND)
      set(_boost_python_suffix_found ${_boost_python_suffix})
      break()
    endif()
  endforeach()

  if(NOT Boost_FOUND)
    message(STATUS "Configuring Python front-end - error (no matching Boost.Python or Boost.NumPy components found)")
  else()
    link_libraries(${Python_LIBRARIES} Boost::python${_boost_python_suffix_found} Boost::numpy${_boost_python_suffix_found} ${TbGAL_LIBRARIES})
    include_directories(${PROJECT_SOURCE_DIR}/python/src ${EIGEN3_INCLUDE_DIRS} ${Python_INCLUDE_DIRS} ${Python_NumPy_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${TbGAL_INCLUDE_DIRS})

    include(GNUInstallDirs)

    file(GLOB module_filenames RELATIVE ${PROJECT_SOURCE_DIR}/python/src ${PROJECT_SOURCE_DIR}/python/src/*.cpp)
    foreach(module_filename ${module_filenames})
      string(REPLACE ".cpp" "" module_name ${module_filename})
      Python_add_library(${module_name} ${PROJECT_SOURCE_DIR}/python/src/${module_filename})
      install(TARGETS ${module_name} LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/tbgal/python/${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/tbgal)
    endforeach()

    install(CODE "file(TOUCH \"${CMAKE_INSTALL_PREFIX}/lib/tbgal/python/${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/tbgal/__init__.py\")")

    message(STATUS "Configuring Python front-end - done")
  endif()
endif()