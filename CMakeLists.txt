cmake_minimum_required(VERSION 3.21)
project(libnurbs VERSION 1.0.0 LANGUAGES CXX)

set(LIBNURBS_SOVERSION 1)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Enable "/utf-8" flag for MSVC compiler
add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

# Add an option to build a shared or static library (default is static)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

find_package(Eigen3 CONFIG REQUIRED)

set(LIBNURBS_INCLUDES
        include/libnurbs/libnurbs.hpp
        include/libnurbs/Curve/Curve.hpp
        include/libnurbs/Core/Typedefs.hpp
        include/libnurbs/Core/KnotVector.hpp
        include/libnurbs/Basis/BSplineBasis.hpp
        include/libnurbs/Algorithm/MathUtils.hpp
        include/libnurbs/Core/Grid.hpp
        include/libnurbs/Surface/Surface.hpp
        include/libnurbs/Geometry/GeomSegment.hpp
        include/libnurbs/Geometry/GeomRect.hpp
        include/libnurbs/Algorithm/KnotRemoval.hpp
        include/libnurbs/Algorithm/DegreeAlgo.hpp
        include/libnurbs/Utils/Serialization.hpp
)

set(LIBNURBS_SOURCES
        src/Core/KnotVector.cpp
        src/Basis/BSplineBasis.cpp
        src/Curve/Curve.cpp
        src/Surface/Surface.cpp
        src/Geometry/GeomSegment.cpp
        src/Geometry/GeomRect.cpp
        src/Algorithm/KnotRemoval.cpp
        src/Algorithm/DegreeAlgo.cpp
        src/Utils/Serialization.cpp
)

# Build library as static or shared based on BUILD_SHARED_LIBS
add_library(libnurbs ${LIBNURBS_SOURCES} ${LIBNURBS_INCLUDES})
add_library(libnurbs::libnurbs ALIAS libnurbs)

# Set the library version and SOVERSION
set_target_properties(libnurbs PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${LIBNURBS_SOVERSION}
        PUBLIC_HEADER "${LIBNURBS_INCLUDES}"
)

# Include directories for the library
target_include_directories(libnurbs
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(libnurbs PUBLIC Eigen3::Eigen)

# Install the library, headers, and cmake configuration files
install(TARGETS libnurbs
        EXPORT libnurbsTargets
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
        PUBLIC_HEADER DESTINATION include/libnurbs
)

# Export target configuration for find_package
install(EXPORT libnurbsTargets
        FILE libnurbsTargets.cmake
        NAMESPACE libnurbs::
        DESTINATION lib/cmake/libnurbs
)

# Generate and install CMake package configuration files
include(CMakePackageConfigHelpers)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfigVersion.cmake
        COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libnurbsConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfig.cmake
        INSTALL_DESTINATION lib/cmake/libnurbs
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libnurbsConfigVersion.cmake
        DESTINATION lib/cmake/libnurbs
)

# Optional build targets
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
