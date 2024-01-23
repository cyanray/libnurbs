cmake_minimum_required(VERSION 3.21)
project(libnurbs_benchmark)
find_package(benchmark CONFIG REQUIRED)

set(libnurbs_Benchmark_SOURCES
        BM_Basis.cpp
)

add_executable(${PROJECT_NAME} ${libnurbs_Benchmark_SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE benchmark::benchmark benchmark::benchmark_main)
target_link_libraries(${PROJECT_NAME} PRIVATE libnurbs)