# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/home/zhipz/work/scripts

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/home/zhipz/work/scripts/build

# Utility rule file for Class.

# Include any custom commands dependencies for this target.
include CMakeFiles/Class.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/Class.dir/progress.make

CMakeFiles/Class: Class.cxx
CMakeFiles/Class: libClass_rdict.pcm
CMakeFiles/Class: libClass.rootmap

Class.cxx: ../include/LinkDef.h
Class.cxx: ../include/Photon_bunches.h
Class.cxx: ../include/events.h
Class.cxx: ../include/Photon_bunches.h
Class.cxx: ../include/events.h
Class.cxx: ../include/LinkDef.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/data/home/zhipz/work/scripts/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating Class.cxx, libClass_rdict.pcm, libClass.rootmap"
	/usr/local/bin/cmake -E env LD_LIBRARY_PATH=/data/home/zhipz/root/lib:/data/home/zhipz/root/lib::/opt/torque-6.1.2/lib:/opt/maui-3.3.1/lib:/data/home/zhipz/root/lib:/data/home/zhipz/hessioxxx/lib:/data/home/zhipz/anaconda3/lib/usr/lib64 /data/home/zhipz/root/bin/rootcling -v2 -f Class.cxx -s /data/home/zhipz/work/scripts/build/libClass.so -rml libClass.so -rmf /data/home/zhipz/work/scripts/build/libClass.rootmap -I/data/home/zhipz/root/include -I/data/home/zhipz/work/scripts /data/home/zhipz/work/scripts/include/Photon_bunches.h /data/home/zhipz/work/scripts/include/events.h /data/home/zhipz/work/scripts/include/LinkDef.h

libClass_rdict.pcm: Class.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libClass_rdict.pcm

libClass.rootmap: Class.cxx
	@$(CMAKE_COMMAND) -E touch_nocreate libClass.rootmap

Class: CMakeFiles/Class
Class: Class.cxx
Class: libClass.rootmap
Class: libClass_rdict.pcm
Class: CMakeFiles/Class.dir/build.make
.PHONY : Class

# Rule to build all files generated by this target.
CMakeFiles/Class.dir/build: Class
.PHONY : CMakeFiles/Class.dir/build

CMakeFiles/Class.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Class.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Class.dir/clean

CMakeFiles/Class.dir/depend:
	cd /data/home/zhipz/work/scripts/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/home/zhipz/work/scripts /data/home/zhipz/work/scripts /data/home/zhipz/work/scripts/build /data/home/zhipz/work/scripts/build /data/home/zhipz/work/scripts/build/CMakeFiles/Class.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Class.dir/depend

