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
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/woosubkim/Documents/workspace/foldseek

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug

# Utility rule file for uninstall.

# Include any custom commands dependencies for this target.
include lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/progress.make

lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall:
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/mmseqs/lib/zstd/build/cmake/lib && /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -P /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/mmseqs/lib/zstd/build/cmake/lib/cmake_uninstall.cmake

uninstall: lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall
uninstall: lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/build.make
.PHONY : uninstall

# Rule to build all files generated by this target.
lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/build: uninstall
.PHONY : lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/build

lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/clean:
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/mmseqs/lib/zstd/build/cmake/lib && $(CMAKE_COMMAND) -P CMakeFiles/uninstall.dir/cmake_clean.cmake
.PHONY : lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/clean

lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/depend:
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/woosubkim/Documents/workspace/foldseek /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/zstd/build/cmake/lib /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/mmseqs/lib/zstd/build/cmake/lib /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/mmseqs/lib/zstd/build/cmake/lib/CMakeFiles/uninstall.dir/depend

