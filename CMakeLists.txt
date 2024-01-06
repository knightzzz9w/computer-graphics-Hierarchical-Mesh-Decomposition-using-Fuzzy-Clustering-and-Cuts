cmake_minimum_required(VERSION 3.0)
project(segmentNet)

# Find PCL
#find_package(PCL 1.8 REQUIRED)

set(CMAKE_BUILD_TYPE "Release")

#set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Eigen is header-only and does not need to be linked
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})


# Include your headers
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR} # or the directory where utils.h is located
    ${CMAKE_CURRENT_SOURCE_DIR}/common
    )

add_library(utils SHARED  common/utils.h common/utils.cpp)

target_link_libraries(utils ${PCL_LIBRARIES} Eigen3::Eigen)

# Add your main application executable
add_executable(${PROJECT_NAME} main.cpp)

set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})
# Link your application with the utils library

target_link_libraries(${PROJECT_NAME} utils)

add_executable(${PROJECT_NAME}_test test.cpp)

set_target_properties(${PROJECT_NAME}_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME}_test utils ${GTEST_LIBRARIES} pthread)