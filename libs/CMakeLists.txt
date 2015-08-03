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


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11" )
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lz -ldl")

set(ARMA_LIB /usr/lib/libarmadillo.so)
set(NLOPT_LIB /usr/local/lib/libnlopt.so)
set(NLOPT_INCLUDE_DIR /usr/local/include)

if (WIN32)
  set(Boost_USE_STATIC_LIBS ON)
endif(WIN32)

add_definitions(${HDF5_DEFINITIONS})

file(GLOB LIBS_FILES "*.?pp" "*.cu")

add_lib_dir(config_IO)
add_lib_dir(data_IO)
add_lib_dir(structure_IO)
add_lib_dir(wavefunctions)
add_lib_dir(potentials)
add_lib_dir(experiments)

include_directories(
  ${PROJECT_SOURCE_DIR}
  config_IO
  data_IO
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

if (CMAKE_BUILD_TYPE MATCHES Release)
        if (MSVC)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/EHsc;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/MD;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/W3;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/nologo;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/O2;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/Zi;)
        else ()
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-O2;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-pipe;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-mtune=generic;)   
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-fstack-protector;) 
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-Wall;) 
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-W;) 
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-D_REENTRANT;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-fPIE;)             
        endif ()
else ()
        if (MSVC)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/EHsc;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/W3;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/nologo;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/Od;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/Zi;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/MDd;)
        else ()
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-O0;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-g;)
        endif ()
endif ()
foreach(cuda_gpu ${CUDA_CAPABILITY_LIST}) 
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-std=c++11;)
list(APPEND CUDA_NVCC_FLAGS "-gencode ${cuda_gpu}")
endforeach(cuda_gpu)

add_library(${PROJECT_NAME} SHARED ${LIBS_FILES})
target_link_libraries( ${PROJECT_NAME} hdfboost ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${FFTW3_LIBS} ${FFTW3F_LIBS} ${FFTW3_OMP_LIBS} ${OPENBABEL_LIB} ${ARMA_LIB} ${NLOPT_LIB}  ${ArrayFire_LIBRARIES} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cufft_LIBRARY} )

if(OPENMP)
  SET_TARGET_PROPERTIES(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${OpenMP_C_FLAGS}" LINK_FLAGS  "${OpenMP_C_FLAGS}")      
endif(OPENMP)
