cmake_minimum_required (VERSION 2.6)
project (hdfboost)
include_directories(
	${CMAKE_SOURCE_DIR}/hdf5boost/lib/
	${ARMA_INCLUDE_DIR}
)
FILE(GLOB antlr_headers2 ${CMAKE_SOURCE_DIR}/hdf5boost/lib/*.hpp)
FILE(GLOB antlr_cfiles ${CMAKE_SOURCE_DIR}/hdf5boost/lib/*.cpp)

add_library(hdfboost SHARED ${antlr_cfiles} ${antlr_headers})
target_link_libraries(hdfboost
${Boost_LIBRARIES}
${HDF5_LIBRARIES}
${ARMA_LIB})
