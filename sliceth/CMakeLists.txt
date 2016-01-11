
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
CMAKE_POLICY(VERSION 2.6)

IF(LUAROCKS_PREFIX)
    MESSAGE(STATUS "Installing Torch through Luarocks")
    STRING(REGEX REPLACE "(.*)lib/luarocks/rocks.*" "\\1" CMAKE_INSTALL_PREFIX  "${LUAROCKS_PREFIX}")
    MESSAGE(STATUS "Prefix inferred from Luarocks: ${CMAKE_INSTALL_PREFIX}")
ENDIF()

FIND_PACKAGE(Torch REQUIRED)

FILE(GLOB luasrc *.lua)
SET(luasrc ${luasrc} test/test.lua)


# INSTALL(FILES ${luasrc}
#   DESTINATION ${ENV_HOME}/torch/install/share/lua/5.1/sliceth)

# didn't work
INSTALL(
  FILES
  ${luasrc}
  DESTINATION "${Torch_INSTALL_LUA_PATH_SUBDIR}/sliceth")

SET(LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../lib)

INSTALL(FILES
  ${LIB_DIR}/exported.hpp
  ${LIB_DIR}/config_IO/configs.hpp
  DESTINATION "${Torch_DIR}/../../../include/sliceth"
)

set(slicelib ../build/debug/bin/libslice.so)
set(hdflib ../build/debug/bin/libhdfboost.so)

if(NOT EXISTS ${slicelib} OR NOT EXISTS ${hdflib})
  MESSAGE(STATUS "Compile the slicepp project in release mode before building sliceth!")
endif(NOT EXISTS ${slicelib} OR NOT EXISTS ${hdflib})

INSTALL(FILES ${slicelib} DESTINATION "${Torch_DIR}/../../../lib/")
INSTALL(FILES ${hdflib} DESTINATION "${Torch_DIR}/../../../lib/")
