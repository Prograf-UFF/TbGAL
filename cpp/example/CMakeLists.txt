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

project(TbGAL-Example)

find_package(TbGAL REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(target_prefix TbGAL_example)

file(GLOB example_filenames RELATIVE ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/*.cpp)
foreach(example_filename ${example_filenames})
  string(REPLACE ".cpp" "" example_name ${example_filename})
  add_executable(${target_prefix}_${example_name} ${PROJECT_SOURCE_DIR}/${example_filename})
  target_include_directories(${target_prefix}_${example_name} PUBLIC ${PROJECT_SOURCE_DIR} ${TbGAL_INCLUDE_DIRS})
  target_link_libraries(${target_prefix}_${example_name})
endforeach()
