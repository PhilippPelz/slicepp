# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /space/projects/slicepp_cuda

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /space/projects/slicepp_cuda

# Include any dependencies generated for this target.
include hdf5boost/lib/CMakeFiles/hdfboost.dir/depend.make

# Include the progress variables for this target.
include hdf5boost/lib/CMakeFiles/hdfboost.dir/progress.make

# Include the compile flags for this target's objects.
include hdf5boost/lib/CMakeFiles/hdfboost.dir/flags.make

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o: hdf5boost/lib/CMakeFiles/hdfboost.dir/flags.make
hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o: hdf5boost/lib/HDFFile.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /space/projects/slicepp_cuda/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o"
	cd /space/projects/slicepp_cuda/hdf5boost/lib && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/hdfboost.dir/HDFFile.cpp.o -c /space/projects/slicepp_cuda/hdf5boost/lib/HDFFile.cpp

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/hdfboost.dir/HDFFile.cpp.i"
	cd /space/projects/slicepp_cuda/hdf5boost/lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /space/projects/slicepp_cuda/hdf5boost/lib/HDFFile.cpp > CMakeFiles/hdfboost.dir/HDFFile.cpp.i

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/hdfboost.dir/HDFFile.cpp.s"
	cd /space/projects/slicepp_cuda/hdf5boost/lib && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /space/projects/slicepp_cuda/hdf5boost/lib/HDFFile.cpp -o CMakeFiles/hdfboost.dir/HDFFile.cpp.s

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.requires:
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.requires

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.provides: hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.requires
	$(MAKE) -f hdf5boost/lib/CMakeFiles/hdfboost.dir/build.make hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.provides.build
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.provides

hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.provides.build: hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o

# Object files for target hdfboost
hdfboost_OBJECTS = \
"CMakeFiles/hdfboost.dir/HDFFile.cpp.o"

# External object files for target hdfboost
hdfboost_EXTERNAL_OBJECTS =

bin/libhdfboost.so: hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o
bin/libhdfboost.so: hdf5boost/lib/CMakeFiles/hdfboost.dir/build.make
bin/libhdfboost.so: /usr/local/lib/libboost_filesystem.so
bin/libhdfboost.so: /usr/local/lib/libboost_system.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libz.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libdl.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libm.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_cpp.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libz.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libdl.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libm.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_hl.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libz.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libdl.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libm.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_hl.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_cpp.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libz.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libdl.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libm.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_cpp.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_hl.so
bin/libhdfboost.so: /usr/local/lib/libhdf5.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libz.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libdl.so
bin/libhdfboost.so: /usr/lib/x86_64-linux-gnu/libm.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_cpp.so
bin/libhdfboost.so: /usr/local/lib/libhdf5_hl.so
bin/libhdfboost.so: hdf5boost/lib/CMakeFiles/hdfboost.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../../bin/libhdfboost.so"
	cd /space/projects/slicepp_cuda/hdf5boost/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/hdfboost.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
hdf5boost/lib/CMakeFiles/hdfboost.dir/build: bin/libhdfboost.so
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/build

hdf5boost/lib/CMakeFiles/hdfboost.dir/requires: hdf5boost/lib/CMakeFiles/hdfboost.dir/HDFFile.cpp.o.requires
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/requires

hdf5boost/lib/CMakeFiles/hdfboost.dir/clean:
	cd /space/projects/slicepp_cuda/hdf5boost/lib && $(CMAKE_COMMAND) -P CMakeFiles/hdfboost.dir/cmake_clean.cmake
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/clean

hdf5boost/lib/CMakeFiles/hdfboost.dir/depend:
	cd /space/projects/slicepp_cuda && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /space/projects/slicepp_cuda /space/projects/slicepp_cuda/hdf5boost/lib /space/projects/slicepp_cuda /space/projects/slicepp_cuda/hdf5boost/lib /space/projects/slicepp_cuda/hdf5boost/lib/CMakeFiles/hdfboost.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : hdf5boost/lib/CMakeFiles/hdfboost.dir/depend

