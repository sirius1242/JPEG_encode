cmake_minimum_required(VERSION 3.18)
project(jpeg_r)

file(GLOB SOURCES ${CMAKE_CURRENT_LIST_DIR}/JPEG.cpp)
add_executable(jpeg_r ${SOURCES})
target_link_libraries(jpeg_r ${MATH_LIBRARIES})