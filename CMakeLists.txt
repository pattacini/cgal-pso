# Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
# Authors: Ugo Pattacini
# CopyPolicy: Released under the terms of the GNU GPL v2.0.

cmake_minimum_required(VERSION 2.8)
project(cgal-pso)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

find_package(YARP)
find_package(ICUBcontrib)
list(APPEND CMAKE_MODULE_PATH ${YARP_MODULE_PATH})
list(APPEND CMAKE_MODULE_PATH ${ICUBCONTRIB_MODULE_PATH})
find_package(CGAL REQUIRED COMPONENTS)
find_package(Boost REQUIRED)

include(ICUBcontribHelpers)
include(ICUBcontribOptions)
icubcontrib_set_default_prefix()

set(sources src/main.cpp src/pso.cpp)
set(headers src/optimizer.h src/pso.h)
source_group("Source Files" FILES ${sources})
source_group("Header Files" FILES ${headers})

include_directories(${YARP_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/src)
include(${CGAL_USE_FILE})
include(CGAL_CreateSingleSourceCGALProgram)

if(MSVC)
  add_definitions(-D_USE_MATH_DEFINES)
endif()

add_executable(${PROJECT_NAME} ${headers} ${sources})
target_link_libraries(${PROJECT_NAME} ${YARP_LIBRARIES})
install(TARGETS ${PROJECT_NAME} DESTINATION bin)
icubcontrib_add_uninstall_target()
