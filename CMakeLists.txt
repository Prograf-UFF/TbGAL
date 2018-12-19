cmake_minimum_required(VERSION 2.8)

find_package(CUDA QUIET REQUIRED)
set(CUDA_HOST_COMPILER /usr/bin/gcc-6)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/TbGAl)
file(WRITE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/__init.py__ "")

set(
    CUDA_NVCC_FLAGS
    ${CUDA_NVCC_FLAGS};
    -O3 --gpu-architecture=compute_35 --gpu-code=sm_35 --expt-relaxed-constexpr -Xcompiler -fPIC -shared -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0/1
    )

include_directories(
    /home/eduardovera/libs/CUSP/
	/home/eduardovera/libs/BOOST/
	/home/eduardovera/.anaconda3/include/python3.7m
    )

link_directories(
	/usr/local/cuda/lib64
	/home/eduardovera/libs/BOOST/lib
	/usr/local/lib
    /usr/lib/x86_64-linux-gnu
	/home/eduardovera/Workspace/TbGAl/build/lib/TbGAl
	/home/eduardovera/.anaconda3/lib
    )

cuda_add_library(metric SHARED src/modules/metric/builder.cu)
SET_TARGET_PROPERTIES(metric PROPERTIES PREFIX "")
target_link_libraries(
    metric
	boost_python37
	cuda
	cudart
	)

cuda_add_library(sparse_tensor SHARED src/modules/sparse_tensor/builder.cu)
SET_TARGET_PROPERTIES(sparse_tensor PROPERTIES PREFIX "")
target_link_libraries(
    sparse_tensor
	boost_python37
	cuda
	cudart
	)

cuda_add_library(multivector SHARED src/modules/multivector/builder.cu)
SET_TARGET_PROPERTIES(multivector PROPERTIES PREFIX "")
target_link_libraries(
    multivector
	boost_python37
	cuda
	cudart
	)

#cuda_add_library(tbgal SHARED src/modules/multivector/multivector.cu)
#SET_TARGET_PROPERTIES(tbgal PROPERTIES PREFIX "")
#target_link_libraries(
#    tbgal
#	boost_python3
#	)

include_directories(
#	src/modules/test
    src/modules/multivector
    src/modules/metric
    src/modules/sparse_tensor
    )

#cuda_add_executable(
#    tests
#    src/modules/test/teste.cu
#    )

#add_dependencies(
#	tests
#	tbgal
#)

#target_link_libraries(
#    tests
#	-lboost_python37
#	-ltbgal
#    -lpython3.7m
#	)
