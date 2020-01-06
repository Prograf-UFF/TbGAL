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

project(TbGAL-Test)

enable_testing()

find_package(TbGAL REQUIRED)
find_package(GATL REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

###############################################################################
# Download and unpack googletest at configure time
configure_file(CMakeLists-googletest.txt.in googletest-download/CMakeLists.txt)

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" . RESULT_VARIABLE result WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} --build . RESULT_VARIABLE result WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src ${CMAKE_CURRENT_BINARY_DIR}/googletest-build EXCLUDE_FROM_ALL)
###############################################################################

set(target_prefix TbGAL_test_)

include_directories(${TbGAL_INCLUDE_DIRS} ${GATL_INCLUDE_DIRS})

file(GLOB test_filenames RELATIVE "${PROJECT_SOURCE_DIR}" ${PROJECT_SOURCE_DIR}/*.cpp)
foreach(test_filename ${test_filenames})
  string(REPLACE ".cpp" "" test_name ${test_filename})
  add_executable(${target_prefix}${test_name} "${PROJECT_SOURCE_DIR}/${test_filename}")
  target_include_directories(${target_prefix}${test_name} PUBLIC ${PROJECT_SOURCE_DIR})
  target_link_libraries(${target_prefix}${test_name} gtest_main ${CMAKE_THREAD_LIBS_INIT})
  add_test(${target_prefix}${test_name} ${target_prefix}${test_name})
endforeach()
