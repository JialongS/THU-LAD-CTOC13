cmake_minimum_required(VERSION 3.14)

project(optimization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(nlopt)


aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/DE DE_SOURCE_FILE)
aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/PSO PSO_SOURCE_FILE)


add_library(optimization ${DE_SOURCE_FILE} ${PSO_SOURCE_FILE} nlopt.hpp nlopt_main.cpp nlopt_main.h rand_safe.h nlopt.h)


target_include_directories(optimization PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/DE
        ${CMAKE_CURRENT_SOURCE_DIR}/PSO
        ${CMAKE_CURRENT_SOURCE_DIR}/ltfo
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${NLOPT_INCLUDE_DIRS})

find_package(OpenMP)
if (OpenMP_FOUND)
    message(STATUS "OpenMP found")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}${OpenMP_CXX_FLAGS}")
endif ()

#find_package(TBB REQUIRED)


target_link_libraries(optimization nlopt OpenMP::OpenMP_CXX)


install(FILES ${CMAKE_INSTALL_PREFIX}/optimization.h DESTINATION ${CMAKE_INSTALL_PREFIX}/optimization)