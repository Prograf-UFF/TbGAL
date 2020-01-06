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
set(VERSION_PATCH 20200104)

project(TbGAL
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
  DESCRIPTION "TbGAL: Tensor-based Geometric Algebra Library"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 NO_MODULE REQUIRED)

file(INSTALL "./cpp/include" DESTINATION "${CMAKE_INSTALL_PREFIX}" FILES_MATCHING PATTERN "*.hpp")

configure_file("./cmake/TbGALConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")

configure_file("./cmake/TbGALConfigVersion.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")
