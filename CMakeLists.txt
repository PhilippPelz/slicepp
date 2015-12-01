cmake_minimum_required(VERSION 2.8)

# Section forces static linking.  Commented because this breaks
#   testing, because Boost makes static linking against Boot Test
#   libraries difficult, and MCS has better things to do than fight it.
#set(CMAKE_USER_MAKE_RULES_OVERRIDE
	#${CMAKE_SOURCE_DIR}/cmake/Modules/c_flag_overrides.cmake)
#set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX
#	${CMAKE_SOURCE_DIR}/cmake/Modules/cxx_flag_overrides.cmake)

project(qstem)
set (VERSION_MAJOR 3)
set (VERSION_MINOR 0)
if (WIN32)

endif(WIN32)
if(WIN32)
	# Squelch Visual studio's warnings about insecure functions - will replace these over time, but must maintain Linux compatibility.
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	set(Boost_USE_STATIC_LIBS ON)
else(WIN32)
	# specify C++0x support for GCC
	ADD_DEFINITIONS(--std=c++11) # Or -std=c++11
endif(WIN32)

set(OPENBABEL_LIB /usr/local/lib/libopenbabel.so)
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lz -ldl -lpthread -fopenmp")
set(CYCLES_CUDA_BINARIES_ARCH sm_52)

# Define the version number for use in projects
add_definitions( -DVERSION="${VERSION_MAJOR}.${VERSION_MINOR}" )
add_definitions( -DBOOST_LOG_DYN_LINK )
add_definitions( -DUSE_THREADS )
add_definitions( -fopenmp )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set (M_LIB "")

set(CLUSTER OFF)

#only for cluster
if(CLUSTER)
	SET(ArrayFire_DIR "$ENV{HOME}/projects/arrayfire/build20")
	SET(BOOST_ROOT "/data/netapp/boost_1_58_0")
	SET(BOOST_INCLUDEDIR "/data/netapp/boost_1_58_0/boost")
	SET(BOOST_LIBRARYDIR "/data/netapp/boost_1_58_0/stage/lib")
	SET(Boost_NO_SYSTEM_PATHS ON)
	set(OPENBABEL_INCLUDE_DIR ~/local/include/openbabel-2.0/)
	set(OPENBABEL_LIB  $ENV{HOME}/local/lib/libopenbabel.so)
	set(ARMA_LIB  $ENV{HOME}/local/lib/libarmadillo.so)
	set(ARMA_INCLUDE_DIR  $ENV{HOME}/local/include)
	set(NLOPT_LIB ~/local/lib/libnlopt.so)
	set(NLOPT_INCLUDE_DIR ~/local/include)
else(CLUSTER)
	set(OPENBABEL_INCLUDE_DIR /usr/local/include/openbabel-2.0/)
	set(OPENBABEL_LIB /usr/local/lib/libopenbabel.so)
	set(ARMA_LIB /usr/lib/libarmadillo.so)
	set(ARMA_INCLUDE_DIR /usr/local/include)
	set(NLOPT_LIB /usr/local/lib/libnlopt.so)
	set(NLOPT_INCLUDE_DIR /usr/local/include)
  SET(BOOST_ROOT "/usr/lib/x86_64-linux-gnu")
endif(CLUSTER)

if (UNIX)
  set (M_LIB "m")
endif(UNIX)
set(HDF5_USE_STATIC_LIBRARIES OFF)

find_package(Boost 1.55 COMPONENTS system filesystem chrono timer log log_setup thread date_time REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(fftw3 REQUIRED)
find_package(CUDA REQUIRED)
find_package(ArrayFire)
include(FindCUDA)
include(CheckComputeCapability)
string(REPLACE " " ";" CUDA_CAPABILITY_LIST ${CUDA_CAPABILITY_INFO})

include_directories(
	${Boost_INCLUDE_DIRS}
	"${FFTW3_INCLUDE_DIRS}"
	${ArrayFire_INCLUDE_DIRS}
	${CUDA_TOOLKIT_INCLUDE}
)


add_subdirectory(libs)
add_subdirectory(stem3)
add_subdirectory(hdf5boost)

#add_subdirectory(gbmaker)
#add_subdirectory(qscRg12)

OPTION( BUILD_TESTS "Set to ON to enable unit test target generation.  Requires Boost Test binary libraries to be installed." OFF )
if (BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)

  # this enables the RUN_TEST target, which runs tests, but doesn't give much info.
  #add_test (NAME TestStem3 COMMAND test_stem3)
  #add_test (NAME TestGBMaker COMMAND test_gbmaker)

  if(WIN32)
    # helpful custom target for verbose CTest output from:
    # http://stackoverflow.com/questions/16158091/cmake-setting-an-environmental-variable-for-ctest-or-otherwise-getting-failed
    if (CMAKE_CONFIGURATION_TYPES)
      add_custom_target(run_tests COMMAND ${CMAKE_CTEST_COMMAND}
	--force-new-ctest-process --output-on-failure
	--build-config "$<CONFIGURATION>")
    else(CMAKE_CONFIGURATION_TYPES)
      add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
	--force-new-ctest-process --output-on-failure)
    endif(CMAKE_CONFIGURATION_TYPES)
  else(WIN32)
    add_custom_target(test COMMAND ${CMAKE_CTEST_COMMAND}
      --force-new-ctest-process --output-on-failure)
  endif(WIN32)
endif(BUILD_TESTS)

OPTION( OPENMP "Set to ON to enable parallel execution using OpenMP" ON )

if(OPENMP)
  find_package(OpenMP REQUIRED)
  if (OPENMP_FOUND)
      set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif(OPENMP)

if (CMAKE_BUILD_TYPE MATCHES Release)
       if (MSVC)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/EHsc;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/MD;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/W3;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/nologo;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/Od;)
	    list(APPEND CUDA_NVCC_FLAGS -Xcompiler;/Zi;)
        else ()
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-O0;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-pipe;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-mtune=generic;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-fstack-protector;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-Wall;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-W;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-D_REENTRANT;)
            list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-fPIE;)
        endif ()
else (CMAKE_BUILD_TYPE MATCHES Release)
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
list(APPEND CUDA_NVCC_FLAGS -Xcompiler;-std=c++11)
list(APPEND CUDA_NVCC_FLAGS "-gencode ${cuda_gpu}")
endforeach(cuda_gpu)

if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Debug")
	set(CMAKE_BUILD_TYPE "Debug")
endif(NOT CMAKE_BUILD_TYPE)
