cmake_minimum_required(VERSION 2.8)

find_package(CUDA QUIET REQUIRED)
set(CUDA_HOST_COMPILER /usr/bin/gcc-6)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/TbGAl)
#set(file "${CMAKE_CURRENT_BINARY_DIR}/file-to-touch")
#file(WRITE filename "message to write"... )
file(WRITE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init.py__ "")

set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 --expt-relaxed-constexpr -Xcompiler -fPIC -shared
    )

include_directories(
    /home/eduardovera/Workspace/CUSP/
	/usr/include/python3.6
    )

link_directories(
	/usr/local/lib
    )

cuda_add_library(metric SHARED src/modules/metric/builder.cu)
SET_TARGET_PROPERTIES(metric PROPERTIES PREFIX "")
target_link_libraries(
    metric
	boost_python3
	)

cuda_add_library(sparse_tensor SHARED src/modules/sparse_tensor/builder.cu)
SET_TARGET_PROPERTIES(sparse_tensor PROPERTIES PREFIX "")
target_link_libraries(
    sparse_tensor
	boost_python3
	)

cuda_add_library(multivector SHARED src/modules/multivector/builder.cu)
SET_TARGET_PROPERTIES(multivector PROPERTIES PREFIX "")
target_link_libraries(
    multivector
	boost_python3
	)
