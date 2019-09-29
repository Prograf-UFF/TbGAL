cmake_minimum_required(VERSION 3.10)

project(TbGAL-Test)

enable_testing()

find_package(TbGAL REQUIRED)
find_package(GATL REQUIRED)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(target_prefix TbGAL_test_)

include_directories(${TbGAL_INCLUDE_DIRS} ${GATL_INCLUDE_DIRS} ${GTEST_INCLUDE_DIRS})

file(GLOB test_filenames RELATIVE "${PROJECT_SOURCE_DIR}" ${PROJECT_SOURCE_DIR}/*.cpp)
foreach(test_filename ${test_filenames})
  string(REPLACE ".cpp" "" test_name ${test_filename})
  add_executable(${target_prefix}${test_name} "${PROJECT_SOURCE_DIR}/${test_filename}")
  target_include_directories(${target_prefix}${test_name} PUBLIC ${PROJECT_SOURCE_DIR})
  target_link_libraries(${target_prefix}${test_name} ${GTEST_BOTH_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  add_test(${target_prefix}${test_name} ${target_prefix}${test_name})
endforeach()
