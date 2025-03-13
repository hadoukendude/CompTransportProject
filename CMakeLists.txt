# CMakeList.txt : CMake project for CompTransportProject, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("CompTransportProject")

# Specify the path to Eigen3
set(Eigen3_DIR "C:\\Users\\Andrew\\vcpkg\\installed\\x64-windows\\share\\eigen3")
find_package (Eigen3 REQUIRED NO_MODULE)

# Add source to this project's executable.
add_executable (CompTransportProject "main.cpp" "main.h" "input.cpp" "input.h" "transport.cpp" "geom.cpp" "geom.h" "transport.h")

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET CompTransportProject PROPERTY CXX_STANDARD 20)
endif()

# TODO: Add tests and install targets if needed.
target_link_libraries (CompTransportProject Eigen3::Eigen)