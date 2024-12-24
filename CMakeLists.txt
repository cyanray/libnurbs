cmake_minimum_required(VERSION 3.28)
project(libnurbs VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable "/utf-8" flag for MSVC compiler
add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

find_package(Eigen3 CONFIG REQUIRED)

option(BUILD_SHARED_LIBS "Build as a shared library" OFF)

add_library(libnurbs)
add_library(libnurbs::libnurbs ALIAS libnurbs)

target_include_directories(libnurbs
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(libnurbs PUBLIC Eigen3::Eigen)

# Add header files to project #
file(GLOB_RECURSE HEADER_FILES
        CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)

target_sources(${PROJECT_NAME}
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS include
        FILES ${HEADER_FILES}
)

# Add source code files to project #
add_subdirectory(src)

set_target_properties(libnurbs PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR})

if (BUILD_SHARED_LIBS)
    set_target_properties(libnurbs PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif ()


install(TARGETS libnurbs
        EXPORT libnurbsTargets
        FILE_SET HEADERS
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfigVersion.cmake
        COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libnurbsConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfig.cmake
        INSTALL_DESTINATION cmake/libnurbs
)


install(EXPORT libnurbsTargets
        FILE libnurbsTargets.cmake
        NAMESPACE libnurbs::
        DESTINATION cmake/libnurbs
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfigVersion.cmake
        DESTINATION cmake/libnurbs
)


option(LIBNURBS_BUILD_UNITTEST "Build Unit-Tests" OFF)
if (LIBNURBS_BUILD_UNITTEST)
    enable_testing()
    add_subdirectory(UnitTest)
endif()

option(LIBNURBS_BUILD_VISUALTEST "Build Visual Tests" OFF)
if (LIBNURBS_BUILD_VISUALTEST)
    add_subdirectory(VisualTest)
endif()

option(LIBNURBS_BUILD_BENCHMARK "Build Benchmark" OFF)
if (LIBNURBS_BUILD_BENCHMARK)
    add_subdirectory(Benchmark)
endif()
