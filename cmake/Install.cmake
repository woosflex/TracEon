include(CMakePackageConfigHelpers)

# Keep bundled codecs private to this installation instead of replacing system libz/liblz4.
get_target_property(TRACEON_ZLIB_TARGET zlibstatic ALIASED_TARGET)
if(NOT TRACEON_ZLIB_TARGET)
    set(TRACEON_ZLIB_TARGET zlibstatic)
endif()
set_target_properties(${TRACEON_ZLIB_TARGET} PROPERTIES
    EXPORT_NAME bundled_zlib POSITION_INDEPENDENT_CODE ON)
set_target_properties(lz4_static PROPERTIES EXPORT_NAME bundled_lz4)
install(TARGETS ${TRACEON_ZLIB_TARGET} lz4_static EXPORT TracEonTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/TracEon COMPONENT TracEon)
install(TARGETS traceon_core traceon_kmer EXPORT TracEonTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TracEon)
install(TARGETS traceon_make_cache
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT TracEon)
set(TRACEON_REMOTE_INSTALLED OFF)
if(TARGET traceon_remote)
    set(TRACEON_REMOTE_INSTALLED ON)
    install(TARGETS traceon_remote EXPORT TracEonTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT TracEon)
    install(TARGETS traceon
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT TracEon)
endif()
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/TracEon COMPONENT TracEon
    FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${ankerl_unordered_dense_SOURCE_DIR}/include/ankerl
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/TracEon COMPONENT TracEon)
install(FILES ${PROJECT_SOURCE_DIR}/LICENSE
    DESTINATION ${CMAKE_INSTALL_DATADIR}/TracEon/licenses COMPONENT TracEon)
install(FILES ${ankerl_unordered_dense_SOURCE_DIR}/LICENSE
    DESTINATION ${CMAKE_INSTALL_DATADIR}/TracEon/licenses RENAME unordered_dense-LICENSE COMPONENT TracEon)
install(FILES ${zlib-ng_SOURCE_DIR}/LICENSE.md
    DESTINATION ${CMAKE_INSTALL_DATADIR}/TracEon/licenses RENAME zlib-ng-LICENSE.md COMPONENT TracEon)
install(FILES ${PROJECT_SOURCE_DIR}/third_party/lz4/lib/LICENSE
    DESTINATION ${CMAKE_INSTALL_DATADIR}/TracEon/licenses RENAME lz4-LICENSE COMPONENT TracEon)

set(TRACEON_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/TracEon)
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/TracEonConfig.cmake.in
    ${PROJECT_BINARY_DIR}/TracEonConfig.cmake INSTALL_DESTINATION ${TRACEON_CONFIG_DIR})
write_basic_package_version_file(${PROJECT_BINARY_DIR}/TracEonConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(EXPORT TracEonTargets NAMESPACE TracEon:: DESTINATION ${TRACEON_CONFIG_DIR} COMPONENT TracEon)
install(FILES ${PROJECT_BINARY_DIR}/TracEonConfig.cmake ${PROJECT_BINARY_DIR}/TracEonConfigVersion.cmake
    DESTINATION ${TRACEON_CONFIG_DIR} COMPONENT TracEon)

if(BUILD_TESTING AND Python3_Interpreter_FOUND)
    add_test(NAME installed_consumer COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/tests/test_install.py ${CMAKE_COMMAND} ${PROJECT_BINARY_DIR}
        ${PROJECT_SOURCE_DIR}/tests/downstream ${CMAKE_CXX_COMPILER} $<CONFIG>)
endif()
