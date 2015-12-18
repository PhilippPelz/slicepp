cmake_minimum_required(VERSION 2.8)

project(slicepp CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/lib
	${CMAKE_SOURCE_DIR}/hdf5boost/lib
	${OPENBABEL_INCLUDE_DIR}
	${HDF5_INCLUDE_DIRS}
	${NLOPT_INCLUDE_DIR}
	${ArrayFire_INCLUDE_DIRS}
	${CUDA_TOOLKIT_INCLUDE}
)

FILE(GLOB SLICE_LIB_HEADERS "${CMAKE_SOURCE_DIR}/lib/*.hpp")
FILE(GLOB SLICE_LIB_HEADERS2 "${CMAKE_SOURCE_DIR}/lib/*/*.hpp")
FILE(GLOB SLICE_LIB_HEADERS3 "${CMAKE_SOURCE_DIR}/hdf5boost/lib/*.hpp")
# FILE(GLOB SLICE_LIB_CU "${CMAKE_SOURCE_DIR}/libs/*/*.cu")
FILE(GLOB SLICEPP_C_FILES "${CMAKE_SOURCE_DIR}/exe/*.cpp")
FILE(GLOB SLICEPP_H_FILES "${CMAKE_SOURCE_DIR}/exe/*.hpp")

add_executable(${PROJECT_NAME}
		${SLICEPP_C_FILES}
		${SLICEPP_H_FILES}
		${SLICE_LIB_HEADERS}
		${SLICE_LIB_HEADERS2}
		${SLICE_LIB_HEADERS3}
		${SLICE_LIB_CU})
#
# # m is libm - math libraries on Unix systems
target_link_libraries(${PROJECT_NAME}
		${Boost_LIBRARIES}
		${HDF5_LIBRARIES}
		${FFTW3F_LIBS}
		${FFTW3_OMP_LIBS}
		${OPENBABEL_LIB}
		${ARMA_LIB}
		${NLOPT_LIB}
	 	${ArrayFire_LIBRARIES}
		afcuda
		${NVVM_LIB}
		hdfboost
		slice)
