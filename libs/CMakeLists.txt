cmake_minimum_required(VERSION 2.8)

project(qstem)

function (add_lib_dir dirname)
	#add_subdirectory(${dirname})
        #exploring this: not ideal to use glob because it defeats some of CMake's ability to detect changes
        #  a better way would be to manually add files in CMakeLists.txt in each folder (still adding them ultimately
        #  to LIBS_FILES
	FILE(GLOB TMP "${dirname}/*.?pp" "${dirname}/*.cu")
	set(LIBS_FILES ${LIBS_FILES} ${TMP} PARENT_SCOPE)
	SOURCE_GROUP(${dirname} FILES ${TMP} PARENT_SCOPE)

endfunction(add_lib_dir)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC " )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lz -ldl")

add_definitions(${HDF5_DEFINITIONS})

file(GLOB LIBS_FILES "*.?pp" "*.cu")

add_lib_dir(config_IO)
add_lib_dir(data_IO)
add_lib_dir(detectors)
add_lib_dir(experiments)
add_lib_dir(potentials)
add_lib_dir(structure_IO)
add_lib_dir(wavefunctions)
add_lib_dir(detectors)

include_directories(
  ${PROJECT_SOURCE_DIR}
  config_IO
  data_IO
  detectors
  structure_IO
  wavefunctions
  potentials
  experiments
  ../hdf5boost/lib
  ${OPENBABEL_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS} 
  ${FFTW3_INCLUDE_DIRS} 
  ${HDF5_INCLUDE_DIR}
  ${NLOPT_INCLUDE_DIR}
  ${ArrayFire_INCLUDE_DIRS}
  ${CUDA_TOOLKIT_INCLUDE}
  )

add_library(${PROJECT_NAME} SHARED ${LIBS_FILES})
target_link_libraries( ${PROJECT_NAME} hdfboost ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW3_LIBS} ${FFTW3F_LIBS} ${FFTW3_OMP_LIBS} ${OPENBABEL_LIB} ${ARMA_LIB} ${NLOPT_LIB}  ${ArrayFire_LIBRARIES} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cufft_LIBRARY} )

if(OPENMP)
  SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS  "${OpenMP_C_FLAGS}")      
endif(OPENMP)
