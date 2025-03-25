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


# Specify the path to Boost
find_package(Boost REQUIRED COMPONENTS math)
# Specify the path to Eigen3  
find_package (Eigen3 REQUIRED)


# Add source to this project's executable.  
add_executable (CompTransportProject "main.cpp" "main.h" "input.cpp" "input.h" "tau.cpp" "geom.cpp" "geom.h" "tau.h" "collprob.cpp" "collprob.h" "matrix.h" "matrix.cpp" "eigenvalue.cpp")  

if (CMAKE_VERSION VERSION_GREATER 3.12)  
 set_property(TARGET CompTransportProject PROPERTY CXX_STANDARD 20)  
endif()  

# TODO: Add tests and install targets if needed.  

target_link_libraries (CompTransportProject Boost::math Eigen3::Eigen)