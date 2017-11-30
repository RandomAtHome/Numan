cmake_minimum_required(VERSION 3.6)
project(Numan)

set(CMAKE_C_STANDARD 11)

set(FIRST first_task.c)
set(SECOND second_task.c)
add_executable(Numan_1 ${FIRST})
add_executable(Numan_2 ${SECOND})
target_link_libraries(Numan_1 m)
target_link_libraries(Numan_2 m)