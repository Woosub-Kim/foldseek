# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.25.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.25.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/woosubkim/Documents/workspace/foldseek

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/woosubkim/Documents/workspace/foldseek/build

# Include any dependencies generated for this target.
include lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/progress.make

# Include the compile flags for this target's objects.
include lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/flags.make

lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o: lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/flags.make
lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o: /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/microtar/microtar.c
lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o: lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && /usr/local/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o -MF CMakeFiles/microtar.dir/microtar.c.o.d -o CMakeFiles/microtar.dir/microtar.c.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/microtar/microtar.c

lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/microtar.dir/microtar.c.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && /usr/local/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/microtar/microtar.c > CMakeFiles/microtar.dir/microtar.c.i

lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/microtar.dir/microtar.c.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && /usr/local/bin/gcc-12 $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/microtar/microtar.c -o CMakeFiles/microtar.dir/microtar.c.s

# Object files for target microtar
microtar_OBJECTS = \
"CMakeFiles/microtar.dir/microtar.c.o"

# External object files for target microtar
microtar_EXTERNAL_OBJECTS =

lib/mmseqs/lib/microtar/libmicrotar.a: lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/microtar.c.o
lib/mmseqs/lib/microtar/libmicrotar.a: lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/build.make
lib/mmseqs/lib/microtar/libmicrotar.a: lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C static library libmicrotar.a"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && $(CMAKE_COMMAND) -P CMakeFiles/microtar.dir/cmake_clean_target.cmake
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/microtar.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/build: lib/mmseqs/lib/microtar/libmicrotar.a
.PHONY : lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/build

lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/clean:
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar && $(CMAKE_COMMAND) -P CMakeFiles/microtar.dir/cmake_clean.cmake
.PHONY : lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/clean

lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/depend:
	cd /Users/woosubkim/Documents/workspace/foldseek/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/woosubkim/Documents/workspace/foldseek /Users/woosubkim/Documents/workspace/foldseek/lib/mmseqs/lib/microtar /Users/woosubkim/Documents/workspace/foldseek/build /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar /Users/woosubkim/Documents/workspace/foldseek/build/lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/mmseqs/lib/microtar/CMakeFiles/microtar.dir/depend

