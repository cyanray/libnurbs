# ports/portfile.cmake

set(SOURCE_PATH "${CMAKE_CURRENT_LIST_DIR}/../..")

if(VCPKG_LIBRARY_LINKAGE STREQUAL "dynamic")
    set(BUILD_SHARED_LIBS_OPTION "-DBUILD_SHARED_LIBS=ON")
else()
    set(BUILD_SHARED_LIBS_OPTION "-DBUILD_SHARED_LIBS=OFF")
endif()

vcpkg_configure_cmake(
        SOURCE_PATH ${SOURCE_PATH}
        PREFER_NINJA
        OPTIONS
        ${BUILD_SHARED_LIBS_OPTION}
)


vcpkg_install_cmake()

vcpkg_copy_pdbs()

vcpkg_fixup_cmake_targets(CONFIG_PATH cmake/libnurbs)

if(EXISTS "${CURRENT_PACKAGES_DIR}/debug/include")
    file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
endif()


vcpkg_install_copyright(FILE_LIST  ${SOURCE_PATH}/LICENSE.txt)
