#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "afcuda" for configuration "Release"
set_property(TARGET afcuda APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(afcuda PROPERTIES
  IMPORTED_LINK_INTERFACE_LIBRARIES_RELEASE ""
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libafcuda.so.3.1.0"
  IMPORTED_SONAME_RELEASE "libafcuda.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS afcuda )
list(APPEND _IMPORT_CHECK_FILES_FOR_afcuda "${_IMPORT_PREFIX}/lib/libafcuda.so.3.1.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
