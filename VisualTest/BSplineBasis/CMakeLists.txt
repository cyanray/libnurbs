cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 23)

add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

find_package(Matplot++ CONFIG REQUIRED)

add_executable(main main.cpp)
target_link_libraries(main PRIVATE libnurbs)
target_link_libraries(main PRIVATE Matplot++::cimg Matplot++::matplot)
