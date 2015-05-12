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

if(WIN32)
	# Squelch Visual studio's warnings about insecure functions - will replace these over time, but must maintain Linux compatibility.
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
else(WIN32)
  # specify C++0x support for GCC
  ADD_DEFINITIONS(-std=c++11) # Or -std=c++11
endif(WIN32)

set(OPENBABEL_LIB /usr/local/lib/libopenbabel.so)
SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lz -ldl -lpthread -fopenmp")

# Define the version number for use in projects
add_definitions( -DVERSION="${VERSION_MAJOR}.${VERSION_MINOR}" )
add_definitions( -DBOOST_LOG_DYN_LINK )
add_definitions( -DUSE_THREADS )
add_definitions( -fopenmp )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set (M_LIB "")
set(OPENBABEL_INCLUDE_DIR /usr/local/include/openbabel-2.0/)
set(OPENBABEL_LIB /usr/local/lib/libopenbabel.so)

if (UNIX)
  set (M_LIB "m")
endif(UNIX)
set(HDF5_USE_STATIC_LIBRARIES OFF)

find_package(Boost COMPONENTS system filesystem chrono timer log log_setup thread date_time REQUIRED)
find_package(HDF5 COMPONENTS C CXX HL REQUIRED)					
find_package(fftw3 REQUIRED)

include_directories(
	${Boost_INCLUDE_DIRS} 
	"${FFTW3_INCLUDE_DIRS}"
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



if (NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Debug")
	set(CMAKE_BUILD_TYPE "Debug")
endif(NOT CMAKE_BUILD_TYPE)