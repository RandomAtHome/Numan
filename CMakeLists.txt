cmake_minimum_required(VERSION 3.6)
project(Numan)

set(CMAKE_C_STANDARD 11)

set(SOURCE_FILES main.c)
add_executable(Numan ${SOURCE_FILES})
target_link_libraries(Numan m)