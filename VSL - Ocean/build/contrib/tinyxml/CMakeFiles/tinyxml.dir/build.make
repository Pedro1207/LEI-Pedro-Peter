# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\CMake\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\CMake\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build"

# Include any dependencies generated for this target.
include contrib/tinyxml/CMakeFiles/tinyxml.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include contrib/tinyxml/CMakeFiles/tinyxml.dir/compiler_depend.make

# Include the progress variables for this target.
include contrib/tinyxml/CMakeFiles/tinyxml.dir/progress.make

# Include the compile flags for this target's objects.
include contrib/tinyxml/CMakeFiles/tinyxml.dir/flags.make

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/flags.make
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/includes_CXX.rsp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj: ../contrib/tinyxml/tinyxml.cpp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj -MF CMakeFiles\tinyxml.dir\tinyxml.cpp.obj.d -o CMakeFiles\tinyxml.dir\tinyxml.cpp.obj -c "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxml.cpp"

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyxml.dir/tinyxml.cpp.i"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxml.cpp" > CMakeFiles\tinyxml.dir\tinyxml.cpp.i

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyxml.dir/tinyxml.cpp.s"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxml.cpp" -o CMakeFiles\tinyxml.dir\tinyxml.cpp.s

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/flags.make
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/includes_CXX.rsp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj: ../contrib/tinyxml/tinyxmlerror.cpp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj -MF CMakeFiles\tinyxml.dir\tinyxmlerror.cpp.obj.d -o CMakeFiles\tinyxml.dir\tinyxmlerror.cpp.obj -c "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlerror.cpp"

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.i"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlerror.cpp" > CMakeFiles\tinyxml.dir\tinyxmlerror.cpp.i

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.s"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlerror.cpp" -o CMakeFiles\tinyxml.dir\tinyxmlerror.cpp.s

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/flags.make
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/includes_CXX.rsp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj: ../contrib/tinyxml/tinyxmlparser.cpp
contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj: contrib/tinyxml/CMakeFiles/tinyxml.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj -MF CMakeFiles\tinyxml.dir\tinyxmlparser.cpp.obj.d -o CMakeFiles\tinyxml.dir\tinyxmlparser.cpp.obj -c "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlparser.cpp"

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.i"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlparser.cpp" > CMakeFiles\tinyxml.dir\tinyxmlparser.cpp.i

contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.s"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml\tinyxmlparser.cpp" -o CMakeFiles\tinyxml.dir\tinyxmlparser.cpp.s

# Object files for target tinyxml
tinyxml_OBJECTS = \
"CMakeFiles/tinyxml.dir/tinyxml.cpp.obj" \
"CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj" \
"CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj"

# External object files for target tinyxml
tinyxml_EXTERNAL_OBJECTS =

contrib/tinyxml/libtinyxml.a: contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxml.cpp.obj
contrib/tinyxml/libtinyxml.a: contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlerror.cpp.obj
contrib/tinyxml/libtinyxml.a: contrib/tinyxml/CMakeFiles/tinyxml.dir/tinyxmlparser.cpp.obj
contrib/tinyxml/libtinyxml.a: contrib/tinyxml/CMakeFiles/tinyxml.dir/build.make
contrib/tinyxml/libtinyxml.a: contrib/tinyxml/CMakeFiles/tinyxml.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libtinyxml.a"
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && $(CMAKE_COMMAND) -P CMakeFiles\tinyxml.dir\cmake_clean_target.cmake
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\tinyxml.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
contrib/tinyxml/CMakeFiles/tinyxml.dir/build: contrib/tinyxml/libtinyxml.a
.PHONY : contrib/tinyxml/CMakeFiles/tinyxml.dir/build

contrib/tinyxml/CMakeFiles/tinyxml.dir/clean:
	cd /d "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" && $(CMAKE_COMMAND) -P CMakeFiles\tinyxml.dir\cmake_clean.cmake
.PHONY : contrib/tinyxml/CMakeFiles/tinyxml.dir/clean

contrib/tinyxml/CMakeFiles/tinyxml.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean" "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\contrib\tinyxml" "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build" "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml" "F:\Main\Uni\LEI - GIT\LEI-Pedro-Peter\VSL - Ocean\build\contrib\tinyxml\CMakeFiles\tinyxml.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : contrib/tinyxml/CMakeFiles/tinyxml.dir/depend
