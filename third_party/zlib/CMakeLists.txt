# This file builds the ZLIB static library from its source files.

cmake_minimum_required(VERSION 3.10)
project(zlib)

# Create a static library named 'zlib_static'
add_library(zlib_static STATIC
        adler32.c
        compress.c
        crc32.c
        deflate.c
        gzclose.c
        gzlib.c
        gzread.c
        gzwrite.c
        infback.c
        inffast.c
        inflate.c
        inftrees.c
        trees.c
        uncompr.c
        zutil.c
)

# Add this directory to the include path so headers can be found
target_include_directories(zlib_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# --- Portability Fix for POSIX systems ---
# Check for unistd.h, which contains declarations for read, write, lseek, etc.
include(CheckIncludeFile)
check_include_file(unistd.h HAVE_UNISTD_H)

# If the header exists, add the HAVE_UNISTD_H definition.
# The zlib source code is written to use this macro to include <unistd.h>
# where available, which resolves the implicit function declaration warnings.
if(HAVE_UNISTD_H)
    target_compile_definitions(zlib_static PUBLIC HAVE_UNISTD_H)
endif()

