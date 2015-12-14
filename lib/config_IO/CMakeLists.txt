cmake_minimum_required(VERSION 2.8)

FILE(GLOB C_FILES "*.cpp")
FILE(GLOB H_FILES "*.hpp")

SOURCE_GROUP(config_IO FILES ${C_FILES} ${H_FILES})

# Add files from this folder to the list of files for the shared library
set (LIBS_C_FILES ${LIBS_C_FILES} ${C_FILES} PARENT_SCOPE})
set (LIBS_C_FILES ${LIBS_H_FILES} ${H_FILES} PARENT_SCOPE})
