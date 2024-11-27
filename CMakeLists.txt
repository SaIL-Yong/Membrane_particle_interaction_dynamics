cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_C_FLAGS "-O3")
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_DEBUG_POSTFIX d)
set(CMAKE_C_FLAGS_DEBUG "-DDEBUG -O0 -g")
set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG -O0 -g")
#message("CMAKE_C_FLAGS_DEBUG is ${CMAKE_C_FLAGS_DEBUG}")
#message("CMAKE_C_FLAGS_RELEASE is ${CMAKE_C_FLAGS_RELEASE}")

#message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
#message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")

#message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

project(MemDynamics VERSION 0.2)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(VersionConfig.h.in VersionConfig.h)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
option(LIBIGL_GLFW "Build target igl::glfw" OFF)
include(libigl)

# Add your project files
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/src/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")
list(REMOVE_ITEM SRC_FILES "${PROJECT_SOURCE_DIR}/src/test.cpp")
add_executable(${PROJECT_NAME} ${INC_FILES} ${SRC_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
if (LIBIGL_GLFW)
  target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
endif()

# add the binary tree to the search path for include files
# so that we will find VersionConfig.h
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
