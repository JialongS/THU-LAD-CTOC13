cmake_minimum_required(VERSION 3.22)
project(inplane)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_C_FLAGS_RELEASE "$ENV{CFLAGS} -O3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -std=c++17 ")
set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)


set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/build/install")


# Add sub-projects
add_subdirectory(dace)
add_subdirectory(optimization)
#find_package(dace 2.0.0 REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CURL REQUIRED libcurl)
include_directories(${CURL_INCLUDE_DIRS})


# Gather source files
aux_source_directory(orbit SOURCE_FILES_ORBIT)
aux_source_directory(Optimization SOURCE_FILE_OPTIMIZATION)

# Add cminpack
file(GLOB MINPACK_FILE MinPack/*.c)
add_library(MINPACK STATIC ${MINPACK_FILE})
target_include_directories(MINPACK PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/MinPack)

# OpenMP
find_package(OpenMP)
if (OpenMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}${OpenMP_CXX_FLAGS}")
endif ()

# Eigen
#find_package(Eigen3 REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIR})


# Add executable
add_executable(inplane ${SOURCE_FILES_ORBIT} ${SOURCE_FILE_OPTIMIZATION} main.cpp
        io.h
        inPlaneOptimization_tym.cpp
        inPlaneOptimization_tym.h
        inPlaneOptimization.h
        inPlaneOptimization.cpp
)


# Set include directories
target_include_directories(inplane PUBLIC orbit MinPack dace ${NLOPT_INCLUDE_DIRS} optimization)


# Link libraries
target_link_libraries(inplane MINPACK dace optimization nlopt)
target_link_libraries(inplane OpenMP::OpenMP_CXX)