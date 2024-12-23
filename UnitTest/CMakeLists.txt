cmake_minimum_required(VERSION 3.21)
project(libnurbs_unit_tests)
find_package(Catch2 CONFIG REQUIRED)

set(libnurbs_UNITTEST_SOURCES
        MathUnitTest.cpp
        KnotVectorUnitTest.cpp
        BSplineBasisUnitTest.cpp
        CurveUnitTest.cpp
        SurfaceUnitTest.cpp
        GeomSegmentUnitTest.cpp
        GeomRectUnitTest.cpp
        GridUnitTest.cpp
)

add_executable(${PROJECT_NAME} ${libnurbs_UNITTEST_SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE Catch2::Catch2WithMain)
target_link_libraries(${PROJECT_NAME} PRIVATE libnurbs)

include(Catch)
catch_discover_tests(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
)