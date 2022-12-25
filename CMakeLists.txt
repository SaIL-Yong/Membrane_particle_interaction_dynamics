cmake_minimum_required(VERSION 3.12)
project(example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)

# Enable the target igl::glfw
igl_include(glfw)


#file(GLOB SRC_FILES *.cpp)
#add_executable(${PROJECT_NAME} ${SRC_FILES})
#target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add your project files
file(GLOB INC_FILES "${PROJECT_SOURCE_DIR}/src/*.h")
file(GLOB SRC_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${INC_FILES} ${SRC_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::core)
if (LIBIGL_GLFW)
  target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw)
endif()

# add the binary tree to the search path for include files
# so that we will find VersionConfig.h
target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
