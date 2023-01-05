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
include lib/foldcomp/CMakeFiles/foldcomp.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/foldcomp/CMakeFiles/foldcomp.dir/progress.make

# Include the compile flags for this target's objects.
include lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make

lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/amino_acid.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o -MF CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o.d -o CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/amino_acid.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/amino_acid.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/amino_acid.cpp > CMakeFiles/foldcomp.dir/src/amino_acid.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/amino_acid.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/amino_acid.cpp -o CMakeFiles/foldcomp.dir/src/amino_acid.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/atom_coordinate.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o -MF CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o.d -o CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/atom_coordinate.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/atom_coordinate.cpp > CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/atom_coordinate.cpp -o CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/database_reader.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o -MF CMakeFiles/foldcomp.dir/src/database_reader.cpp.o.d -o CMakeFiles/foldcomp.dir/src/database_reader.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/database_reader.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/database_reader.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/database_reader.cpp > CMakeFiles/foldcomp.dir/src/database_reader.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/database_reader.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/database_reader.cpp -o CMakeFiles/foldcomp.dir/src/database_reader.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/discretizer.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o -MF CMakeFiles/foldcomp.dir/src/discretizer.cpp.o.d -o CMakeFiles/foldcomp.dir/src/discretizer.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/discretizer.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/discretizer.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/discretizer.cpp > CMakeFiles/foldcomp.dir/src/discretizer.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/discretizer.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/discretizer.cpp -o CMakeFiles/foldcomp.dir/src/discretizer.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/foldcomp.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o -MF CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o.d -o CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/foldcomp.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/foldcomp.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/foldcomp.cpp > CMakeFiles/foldcomp.dir/src/foldcomp.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/foldcomp.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/foldcomp.cpp -o CMakeFiles/foldcomp.dir/src/foldcomp.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/nerf.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o -MF CMakeFiles/foldcomp.dir/src/nerf.cpp.o.d -o CMakeFiles/foldcomp.dir/src/nerf.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/nerf.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/nerf.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/nerf.cpp > CMakeFiles/foldcomp.dir/src/nerf.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/nerf.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/nerf.cpp -o CMakeFiles/foldcomp.dir/src/nerf.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/sidechain.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o -MF CMakeFiles/foldcomp.dir/src/sidechain.cpp.o.d -o CMakeFiles/foldcomp.dir/src/sidechain.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/sidechain.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/sidechain.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/sidechain.cpp > CMakeFiles/foldcomp.dir/src/sidechain.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/sidechain.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/sidechain.cpp -o CMakeFiles/foldcomp.dir/src/sidechain.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/torsion_angle.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o -MF CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o.d -o CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/torsion_angle.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/torsion_angle.cpp > CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/torsion_angle.cpp -o CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.s

lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/flags.make
lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o: /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/utility.cpp
lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o: lib/foldcomp/CMakeFiles/foldcomp.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o -MF CMakeFiles/foldcomp.dir/src/utility.cpp.o.d -o CMakeFiles/foldcomp.dir/src/utility.cpp.o -c /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/utility.cpp

lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/foldcomp.dir/src/utility.cpp.i"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/utility.cpp > CMakeFiles/foldcomp.dir/src/utility.cpp.i

lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/foldcomp.dir/src/utility.cpp.s"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && /usr/local/bin/g++-12 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp/src/utility.cpp -o CMakeFiles/foldcomp.dir/src/utility.cpp.s

# Object files for target foldcomp
foldcomp_OBJECTS = \
"CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o" \
"CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o" \
"CMakeFiles/foldcomp.dir/src/database_reader.cpp.o" \
"CMakeFiles/foldcomp.dir/src/discretizer.cpp.o" \
"CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o" \
"CMakeFiles/foldcomp.dir/src/nerf.cpp.o" \
"CMakeFiles/foldcomp.dir/src/sidechain.cpp.o" \
"CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o" \
"CMakeFiles/foldcomp.dir/src/utility.cpp.o"

# External object files for target foldcomp
foldcomp_EXTERNAL_OBJECTS =

lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/amino_acid.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/atom_coordinate.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/database_reader.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/discretizer.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/foldcomp.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/nerf.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/sidechain.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/torsion_angle.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/src/utility.cpp.o
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/build.make
lib/foldcomp/libfoldcomp.a: lib/foldcomp/CMakeFiles/foldcomp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/woosubkim/Documents/workspace/foldseek/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX static library libfoldcomp.a"
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && $(CMAKE_COMMAND) -P CMakeFiles/foldcomp.dir/cmake_clean_target.cmake
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/foldcomp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/foldcomp/CMakeFiles/foldcomp.dir/build: lib/foldcomp/libfoldcomp.a
.PHONY : lib/foldcomp/CMakeFiles/foldcomp.dir/build

lib/foldcomp/CMakeFiles/foldcomp.dir/clean:
	cd /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp && $(CMAKE_COMMAND) -P CMakeFiles/foldcomp.dir/cmake_clean.cmake
.PHONY : lib/foldcomp/CMakeFiles/foldcomp.dir/clean

lib/foldcomp/CMakeFiles/foldcomp.dir/depend:
	cd /Users/woosubkim/Documents/workspace/foldseek/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/woosubkim/Documents/workspace/foldseek /Users/woosubkim/Documents/workspace/foldseek/lib/foldcomp /Users/woosubkim/Documents/workspace/foldseek/build /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp /Users/woosubkim/Documents/workspace/foldseek/build/lib/foldcomp/CMakeFiles/foldcomp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/foldcomp/CMakeFiles/foldcomp.dir/depend

