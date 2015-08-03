cmake_minimum_required(VERSION 2.8)

project(stem3)		

include_directories(
	${CMAKE_SOURCE_DIR}/libs
	${CMAKE_SOURCE_DIR}/hdf5boost/lib
	${OPENBABEL_INCLUDE_DIR}
	${Boost_INCLUDE_DIRS} 
  	${FFTW3_INCLUDE_DIRS} 
  	${HDF5_INCLUDE_DIR}
  	${NLOPT_INCLUDE_DIR}
  	${ArrayFire_INCLUDE_DIRS}
	${CUDA_TOOLKIT_INCLUDE}
)


FILE(GLOB QSTEM_LIB_HEADERS "${CMAKE_SOURCE_DIR}/libs/*.hpp")
FILE(GLOB QSTEM_LIB_HEADERS2 "${CMAKE_SOURCE_DIR}/libs/*/*.hpp")
FILE(GLOB QSTEM_LIB_HEADERS3 "${CMAKE_SOURCE_DIR}/hdf5boost/lib/*.hpp")
FILE(GLOB QSTEM_LIB_CU "${CMAKE_SOURCE_DIR}/libs/*/*.cu")
FILE(GLOB STEM3_C_FILES "${CMAKE_SOURCE_DIR}/stem3/*.cpp")
FILE(GLOB STEM3_H_FILES "${CMAKE_SOURCE_DIR}/stem3/*.hpp")


cuda_add_executable(stem3 ${STEM3_C_FILES} ${STEM3_H_FILES} ${QSTEM_LIB_HEADERS} ${QSTEM_LIB_HEADERS2} ${QSTEM_LIB_HEADERS3} ${QSTEM_LIB_CU})

# m is libm - math libraries on Unix systems
target_link_libraries(stem3 qstem hdfboost ${FFTW3_LIBS} ${FFTW3_OMP_LIBS} ${M_LIB} ${Boost_LIBRARIES} ${OPENBABEL_LIB} ${ArrayFire_LIBRARIES} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cufft_LIBRARY})


