cmake_minimum_required(VERSION 3.10)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 20190905)

project(TbGAL
  VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}
  DESCRIPTION "TbGAL: Tensor-based Geometric Algebra Library"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(INSTALL "./cpp/include" DESTINATION "${CMAKE_INSTALL_PREFIX}" FILES_MATCHING PATTERN "*.hpp")

configure_file("./cmake/TbGALConfig.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfig.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")

configure_file("./cmake/TbGALConfigVersion.cmake.in" "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" @ONLY NEWLINE_STYLE UNIX)
file(INSTALL "${CMAKE_CURRENT_BINARY_DIR}/TbGALConfigVersion.cmake" DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/tbgal")
