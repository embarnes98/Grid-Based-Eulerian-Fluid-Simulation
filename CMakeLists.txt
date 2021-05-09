cmake_minimum_required(VERSION 3.12)
#-------------------------------------------------------------------------------------------
# I'm going to use vcpk in most cases for our install of 3rd party libs
# this is going to check the environment variable for CMAKE_TOOLCHAIN_FILE and this must point to where
# vcpkg.cmake is in the University this is set in your .bash_profile to
# export CMAKE_TOOLCHAIN_FILE=/public/devel/2020/vcpkg/scripts/buildsystems/vcpkg.cmake
#-------------------------------------------------------------------------------------------
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{CMAKE_TOOLCHAIN_FILE})
   set(CMAKE_TOOLCHAIN_FILE $ENV{CMAKE_TOOLCHAIN_FILE})
endif()

# Name of the project
project(fluid_sim_sfml_build)
# This is the name of the Exe change this and it will change everywhere
set(TargetName Fluid-Sim-SFML)
# find Qt libs
find_package(SFML COMPONENTS system window graphics REQUIRED)
	# use C++ 17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
# Set the name of the executable we want to build
add_executable(${TargetName})
include_directories(include)
# Add NGL include path
target_sources(${TargetName} PRIVATE 
${PROJECT_SOURCE_DIR}/src/main.cpp
${PROJECT_SOURCE_DIR}/src/Fluid-Sim-SFML.cpp
)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	add_compile_definitions(_USE_MATH_DEFINES)
	add_compile_definitions(NOMINMAX)
	target_compile_definitions(${TargetName} PRIVATE /wd4251)
endif()

# add exe and link libs that must be after the other defines
target_link_libraries(${TargetName} PRIVATE FLAC OpenAL OpenGL Vorbis sfml-graphics)

find_package(GTest CONFIG REQUIRED)
include(GoogleTest) 
enable_testing()
add_executable(Fluid-Sim-SFML-Tests)
target_sources(Fluid-Sim-SFML-Tests PRIVATE tests/Fluid-Sim-SFML-Tests.cpp src/Fluid-Sim-SFML.cpp include/Fluid-Sim-SFML.h)
target_link_libraries(Fluid-Sim-SFML-Tests PRIVATE GTest::gtest GTest::gtest_main )
target_link_libraries(Fluid-Sim-SFML-Tests PRIVATE FLAC OpenAL OpenGL Vorbis sfml-graphics)

gtest_discover_tests(Fluid-Sim-SFML-Tests)