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

# Include any dependencies generated for this target.
include lib/tmalign/CMakeFiles/tmalign.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/tmalign/CMakeFiles/tmalign.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/tmalign/CMakeFiles/tmalign.dir/progress.make

# Include the compile flags for this target's objects.
include lib/tmalign/CMakeFiles/tmalign.dir/flags.make

lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o: lib/tmalign/CMakeFiles/tmalign.dir/flags.make
lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o: ../lib/tmalign/affineneedlemanwunsch.cpp
lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o: lib/tmalign/CMakeFiles/tmalign.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o -MF CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o.d -o CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/affineneedlemanwunsch.cpp

lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/affineneedlemanwunsch.cpp > CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.i

lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/affineneedlemanwunsch.cpp -o CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.s

lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o: lib/tmalign/CMakeFiles/tmalign.dir/flags.make
lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o: ../lib/tmalign/TMalign.cpp
lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o: lib/tmalign/CMakeFiles/tmalign.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o -MF CMakeFiles/tmalign.dir/TMalign.cpp.o.d -o CMakeFiles/tmalign.dir/TMalign.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/TMalign.cpp

lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tmalign.dir/TMalign.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/TMalign.cpp > CMakeFiles/tmalign.dir/TMalign.cpp.i

lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tmalign.dir/TMalign.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign/TMalign.cpp -o CMakeFiles/tmalign.dir/TMalign.cpp.s

# Object files for target tmalign
tmalign_OBJECTS = \
"CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o" \
"CMakeFiles/tmalign.dir/TMalign.cpp.o"

# External object files for target tmalign
tmalign_EXTERNAL_OBJECTS =

lib/tmalign/libtmalign.a: lib/tmalign/CMakeFiles/tmalign.dir/affineneedlemanwunsch.cpp.o
lib/tmalign/libtmalign.a: lib/tmalign/CMakeFiles/tmalign.dir/TMalign.cpp.o
lib/tmalign/libtmalign.a: lib/tmalign/CMakeFiles/tmalign.dir/build.make
lib/tmalign/libtmalign.a: lib/tmalign/CMakeFiles/tmalign.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library libtmalign.a"
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && $(CMAKE_COMMAND) -P CMakeFiles/tmalign.dir/cmake_clean_target.cmake
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tmalign.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/tmalign/CMakeFiles/tmalign.dir/build: lib/tmalign/libtmalign.a
.PHONY : lib/tmalign/CMakeFiles/tmalign.dir/build

lib/tmalign/CMakeFiles/tmalign.dir/clean:
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign && $(CMAKE_COMMAND) -P CMakeFiles/tmalign.dir/cmake_clean.cmake
.PHONY : lib/tmalign/CMakeFiles/tmalign.dir/clean

lib/tmalign/CMakeFiles/tmalign.dir/depend:
	cd /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/woosubkim/Documents/workspace/foldseek /Users/woosubkim/Documents/workspace/foldseek/lib/tmalign /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign /Users/woosubkim/Documents/workspace/foldseek/cmake-build-debug/lib/tmalign/CMakeFiles/tmalign.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/tmalign/CMakeFiles/tmalign.dir/depend

