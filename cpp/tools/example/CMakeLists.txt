cmake_minimum_required(VERSION 3.10)

project(TbGAL-Example)

find_package(TbGAL REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(target_prefix TbGAL_example_)

set(example_groups
  basic
)

include_directories(${TbGAL_INCLUDE_DIRS})

foreach(example_group ${example_groups})
  file(GLOB example_filenames RELATIVE "${PROJECT_SOURCE_DIR}/${example_group}" ${PROJECT_SOURCE_DIR}/${example_group}/*.cpp)
  foreach(example_filename ${example_filenames})
    string(REPLACE ".cpp" "" example_name ${example_filename})
    add_executable(${target_prefix}${example_group}_${example_name} "${PROJECT_SOURCE_DIR}/${example_group}/${example_filename}")
    target_include_directories(${target_prefix}${example_group}_${example_name} PUBLIC ${PROJECT_SOURCE_DIR}/${example_group})
    target_link_libraries(${target_prefix}${example_group}_${example_name})
  endforeach()
endforeach()
