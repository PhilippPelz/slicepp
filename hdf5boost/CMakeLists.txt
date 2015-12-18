cmake_minimum_required (VERSION 2.8)
project (hdf5boost)

set( HDF5_USE_STATIC_LIBRARIES OFF)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++11" )

ADD_DEFINITIONS(-std=c++11)

add_subdirectory(lib)

include_directories(
${CMAKE_SOURCE_DIR}
${CMAKE_SOURCE_DIR}/hdf5boost/lib
${HDF5_INCLUDE_DIRS}
${Boost_INCLUDE_DIRS}
)

add_definitions(${HDF5_DEFINITIONS})
SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lz -ldl")

add_executable(hdfboosttest main.cpp)
target_link_libraries(hdfboosttest hdfboost
${Boost_LIBRARIES}
${HDF5_LIBRARIES}
)
