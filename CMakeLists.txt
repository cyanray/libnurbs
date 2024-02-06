cmake_minimum_required(VERSION 3.21)
project(libnurbs)
set(CMAKE_CXX_STANDARD 23)

add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

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

)

set(LIBNURBS_SOURCES
        src/Core/KnotVector.cpp
        src/Basis/BSplineBasis.cpp
        src/Curve/Curve.cpp
        src/Surface/Surface.cpp
        src/Geometry/GeomSegment.cpp
        src/Geometry/GeomRect.cpp

)

add_library(libnurbs STATIC ${LIBNURBS_INCLUDES} ${LIBNURBS_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC include)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

install(TARGETS ${PROJECT_NAME})

option(LIBNURBS_BUILD_UNITTEST "Build Unit-Tests" OFF)
if (LIBNURBS_BUILD_UNITTEST)
    enable_testing()
    add_subdirectory(UnitTest)
endif ()

option(LIBNURBS_BUILD_VISUALTEST "Build Visual Tests" OFF)
if (LIBNURBS_BUILD_VISUALTEST)
    add_subdirectory(VisualTest)
endif ()

option(LIBNURBS_BUILD_BENCHMARK "Build Benchmark" OFF)
if (LIBNURBS_BUILD_BENCHMARK)
    add_subdirectory(Benchmark)
endif ()