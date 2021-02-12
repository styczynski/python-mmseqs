cmake_minimum_required(VERSION 2.8.12)
project(mmseqs_native)

set(CMAKE_VERBOSE_MAKEFILE ON)

if (FORCE_X86)
    message("-- Force compilation for x86")
    set (CMAKE_C_FLAGS "-m32")
endif ()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.15)

set(Python_ADDITIONAL_VERSIONS
    3.8 3.8.0 3.8.1 3.8.2 3.8.3 3.8.4 3.8.5 3.8.6 3.8.7
    3.9 3.9.0 3.9.1
    3.7 3.7.0 3.7.1 3.7.2 3.7.3 3.7.4 3.7.5 3.7.6 3.7.7 3.7.8 3.7.9
    3.6 3.6.0 3.6.1 3.6.2 3.6.3 3.6.4 3.6.5 3.6.6 3.6.7 3.6.8 3.6.9 3.6.10 3.6.11 3.6.12
    3.5 3.5.0 3.5.1 3.5.2 3.5.3 3.5.4 3.5.5 3.5.6 3.5.7 3.5.8 3.5.9 3.5.10
    3.4 3.4.0 3.4.1 3.4.2 3.4.3 3.4.4 3.4.5 3.4.6 3.4.7 3.4.8 3.4.9 3.4.10 3.3.0
    3.3 3.3.1 3.3.2 3.3.3 3.3.4 3.3.5 3.3.6 3.3.7
    3.2 3.2.0 3.2.1 3.2.2 3.2.3 3.2.4 3.2.5 3.2.6
    3.1 3.1.0 3.1.1 3.1.2 3.1.3 3.1.4 3.1.5
    3.0 3.0.0 3.0.1)

set(CMAKE_BUILD_TYPE Debug)
if (MSVC)
    add_compile_options(/DEBUG:FULL)
else ()
    add_compile_options(-g)
endif ()

SET(SOURCE_DIR "src/mmseqs_native")

message("-- Source Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("-- Project Directory: ${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/lib/cppfs/cmake")

include(AppendTargetProperty)
include(GenerateExportHeader)
include(Custom)
include(HealthCheck)

set(FRAMEWORK_ONLY 0 CACHE BOOL "Framework mode (don't create mmseqs executable)")
set(HAVE_SANITIZER 0 CACHE BOOL "Have sanitizers")
set(INSTALL_UTIL 1 CACHE BOOL "Install utility scripts")
set(VERSION_OVERRIDE "" CACHE STRING "Override version string in help and usage messages")
set(DISABLE_IPS4O 0 CACHE BOOL "Disabling IPS4O sorting library requiring 128-bit compare exchange operations")
set(HAVE_AVX2 0 CACHE BOOL "Have CPU with AVX2")
set(HAVE_SSE4_1 0 CACHE BOOL "Have CPU with SSE4.1")
set(HAVE_SSE2 0 CACHE BOOL "Have CPU with SSE2")
set(HAVE_POWER9 0 CACHE BOOL "Have POWER9 CPU")
set(HAVE_POWER8 0 CACHE BOOL "Have POWER8 CPU")
set(HAVE_ARM8 0 CACHE BOOL "Have ARMv8 CPU")
set(NATIVE_ARCH 1 CACHE BOOL "Assume native architecture for SIMD. Use one of the HAVE_* options or set CMAKE_CXX_FLAGS to the appropriate flags if you disable this.")

set(Python3_FIND_STRATEGY LOCATION)
if(defined $ENV{pythonVersion})
    message("-- Detected pythonVersion env variable. Looking for python ${pythonVersion}")
    set(Python3_ROOT_DIR ${pythonVersion})
    find_package(Python3 ${pythonVersion} REQUIRED COMPONENTS NumPy)
else ()
    find_package(Python3 REQUIRED COMPONENTS NumPy)
    if (Python3_NumPy_FOUND)
        message("-- Detected NumPy installation")
    else ()
        message("-- No NumPy")
    endif ()
endif ()

if (HAVE_SANITIZER)
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif ()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

# find compiler
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message("-- Compiler is clang(++)")
    set(CMAKE_COMPILER_IS_CLANG 1)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("-- Compiler is GNU ")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message("-- Compiler is icc(++)")
    set(CMAKE_COMPILER_IS_ICC 1)
endif ()

# see https://wiki.debian.org/ArchitectureSpecificsMemo for char signedness
set(MMSEQS_CXX_FLAGS "-fsigned-char")

# SIMD instruction sets support
set(MMSEQS_ARCH "")
if (HAVE_AVX2)
    if (CMAKE_COMPILER_IS_CLANG)
        set(MMSEQS_ARCH "${MMSEQS_ARCH} -mavx2 -mcx16")
    else ()
        set(MMSEQS_ARCH "${MMSEQS_ARCH} -mavx2 -mcx16 -Wa,-q")
    endif ()
    set(X64 1)
elseif (HAVE_SSE4_1)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -msse4.1 -mcx16")
    set(X64 1)
elseif (HAVE_SSE2)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -msse2")
    set(DISABLE_IPS4O 1)
    set(X64 1)
elseif (HAVE_POWER9)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -mcpu=power9 -mvsx")
    set(PPC64 1)
elseif (HAVE_POWER8)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -mcpu=power8 -mvsx")
    set(PPC64 1)
elseif (HAVE_ARM8)
    set(MMSEQS_ARCH "${MMSEQS_ARCH} -march=armv8-a+simd")
    set(ARM 1)
endif ()

if (NATIVE_ARCH AND (MMSEQS_ARCH STREQUAL ""))
    if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*|aarch64.*|AARCH64.*)")
        set(ARM 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
        set(PPC64 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "amd64|AMD64")
        set(X64 1)
    elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "x86|X86")
        set(X86 1)
    else ()
        message(WARNING "CPU without native SIMD instructions. Performance will be bad.")
    endif ()
    if (EMSCRIPTEN)
        set(MMSEQS_ARCH "-msimd128 -s WASM=1 -s ASSERTIONS=1")
    elseif (X86 OR X64)
        include(CheckSSEFeatures)
        if (NOT HAVE_SSE4_1_EXTENSIONS)
            if (NOT HAVE_SSE2_EXTENSIONS)
                message(WARNING "At least SSE2 is needed to compile")
            endif ()
            message(WARNING "At least SSE4.1 is needed for best performance")
        endif ()
        if (PPC64)
            set(MMSEQS_ARCH "-mcpu=native")
        else ()
            # clang has a problem with march=native on travis
            if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.0.0")
                set(MMSEQS_ARCH "${SSE_FLAGS}")
            else()
                set(MMSEQS_ARCH "-march=native")
            endif()
        endif ()
    else ()
        if (PPC64)
            set(MMSEQS_ARCH "-mcpu=native")
        else ()
            set(MMSEQS_ARCH "-march=native")
        endif ()
    endif ()
endif ()

if (NOT (MMSEQS_ARCH STREQUAL ""))
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} ${MMSEQS_ARCH}")
endif ()

if (CYGWIN OR ARM OR PPC64)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D_GNU_SOURCE=1")
endif ()

if (CMAKE_COMPILER_IS_ICC)
    # default -fp-model results in inconsistent results in profile search
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -fp-model precise")
endif ()

# Apple specific features
if (APPLE)
    # macOS SDK started using _Atomic (in ucred.h) which g++ does not support
    # __APPLE_API_STRICT_CONFORMANCE makes sysctl.h not include apis like ucred.h
    # See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89864
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D__APPLE_API_STRICT_CONFORMANCE")
endif ()

if (CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D_WITH_GETLINE")
endif ()

if (CYGWIN)
    # default cygwin allocator (dlmalloc) locks on every allocation and destroys MT performance
    add_subdirectory(lib/nedmalloc)
    # IPS4O seems to deadlock on cygwin
    set(DISABLE_IPS4O 1)
endif()

if ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "5.0.0")
        OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9.1.0"))
    # clang before v5 throws compile errors on ips4o
    set(DISABLE_IPS4O 1)
endif ()

if (PPC64)
    # FIXME: investigate why on ppc the regression seems to fail randomly
    set(DISABLE_IPS4O 1)
endif ()

set(MMSEQS_C_FLAGS "${MMSEQS_CXX_FLAGS}")
set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -std=c++17")
# Compiler-specific features
if (CMAKE_COMPILER_IS_CLANG AND (NOT EMSCRIPTEN))
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++17")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -stdlib=libc++")
endif ()


# zstd
# We use ZSTD_findDecompressedSize which is only available with ZSTD_STATIC_LINKING_ONLY
# Thus we cannot use a system provided libzstd
set(ZSTD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd")
set(CMAKE_INSTALL_LIBDIR bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd/build_utils/cmake/CMakeModules")
option(ZSTD_LEGACY_SUPPORT "LEGACY SUPPORT" OFF)
option(ZSTD_BUILD_STATIC "BUILD STATIC LIBRARIES" ON)
option(ZSTD_BUILD_SHARED "BUILD SHARED LIBRARIES" OFF)
option(ZSTD_MULTITHREAD_SUPPORT "MULTITHREADING SUPPORT" OFF)
option(ZSTD_BUILD_PROGRAMS "BUILD PROGRAMS" OFF)
option(ZSTD_BUILD_CONTRIB "BUILD CONTRIB" OFF)
option(ZSTD_BUILD_TESTS "BUILD TESTS" OFF)
include_directories(lib/zstd/lib)
add_subdirectory(lib/zstd/build_utils/cmake/lib EXCLUDE_FROM_ALL)
set_target_properties(libzstd_static PROPERTIES COMPILE_FLAGS "${MMSEQS_C_FLAGS}" LINK_FLAGS "${MMSEQS_C_FLAGS}")

# global utils
include_directories(src/mmseqs_native/src)

# cppfs
add_subdirectory(lib/cppfs)

# spdlog
add_subdirectory(lib/spdlog)

# tinyexpr
include_directories(lib/tinyexpr)
add_subdirectory(lib/tinyexpr EXCLUDE_FROM_ALL)

# microtar
include_directories(lib/microtar)
add_subdirectory(lib/microtar)

# simde
include_directories(lib/simde)

include_directories(lib)
include_directories(lib)
include_directories(lib/gzstream)
include_directories(lib/alp)
include_directories(lib/cacode)
include_directories(lib/ksw2)
include_directories(lib/xxhash)
if (NOT DISABLE_IPS4O)
    include_directories(lib/ips4o)
endif ()

add_subdirectory(lib/cacode)
add_subdirectory(lib/alp)
add_subdirectory(lib/ksw2)
add_subdirectory(data)

add_subdirectory(src/mmseqs_native/src/mmseqs/alignment)
add_subdirectory(src/mmseqs_native/src/mmseqs/clustering)
add_subdirectory(src/mmseqs_native/src/mmseqs/commons)
add_subdirectory(src/mmseqs_native/src/mmseqs/linclust)
add_subdirectory(src/mmseqs_native/src/mmseqs/multihit)
add_subdirectory(src/mmseqs_native/src/mmseqs/prefiltering)
add_subdirectory(src/mmseqs_native/src/mmseqs/taxonomy)
add_subdirectory(src/mmseqs_native/src/mmseqs/util)
add_subdirectory(src/mmseqs_native/src/mmseqs/api)
add_subdirectory(src/mmseqs_native/src/mmseqs/workflow)

# if (NOT FRAMEWORK_ONLY AND INSTALL_UTIL)
#     add_subdirectory(src/mmseqs_native/src/util)
# endif ()

set(MMSEQ_FRAMEWORK_SOURCES
    ${alignment_header_files}
    ${alignment_source_files}
    ${clustering_header_files}
    ${clustering_source_files}
    ${api_header_files}
    ${api_source_files}
    ${commons_header_files}
    ${commons_source_files}
    ${prefiltering_header_files}
    ${prefiltering_source_files}
    ${multihit_header_files}
    ${multihit_source_files}
    ${taxonomy_header_files}
    ${taxonomy_source_files}
    ${linclust_source_files}
    ${util_header_files}
    ${util_source_files}
    ${workflow_source_files}
)
set(MMSEQ_TARGET_SOURCES
    $<TARGET_OBJECTS:alp>
    $<TARGET_OBJECTS:ksw2>
    $<TARGET_OBJECTS:cacode>
)
list(TRANSFORM MMSEQ_FRAMEWORK_SOURCES PREPEND src/mmseqs_native/src/)
set(BINDINGS_COMMON_SRCS
    ${MMSEQ_FRAMEWORK_SOURCES}
    ${MMSEQ_TARGET_SOURCES}
)

# Generate python module

set(MMSEQ_TARGET_INCLUDES
    ${Python3_NumPy_INCLUDE_DIRS}
    ${CMAKE_BINARY_DIR}/generated
    ${PROJECT_BINARY_DIR}/generated
    "src/mmseqs_native/src)"
    "lib/toml11"
    "lib/spdlog/include"
    "lib/cppfs/include"
    "src"
    .
)

set(MMSEQ_TARGET_DEPENDENCIES
    "generated"
)

set(MMSEQ_TARGET_LINKS
    "spdlog"
    "cppfs"
    "Python3::NumPy"
)

 # this will probably break if you omit the "s
string(REPLACE " " ";" MMSEQS_CXX_FLAGS_LIST "${MMSEQS_CXX_FLAGS}")

add_subdirectory(lib/pybind11)
pybind11_add_module(mmseqs_native ${SOURCES} "${SOURCE_DIR}/bindings.cpp" ${BINDINGS_COMMON_SRCS})
#set_target_properties(mmseqs_native PROPERTIES LINK_FLAGS "-Wl,-undefined,error")
target_include_directories(mmseqs_native PUBLIC ${MMSEQ_TARGET_INCLUDES})
add_dependencies(mmseqs_native ${MMSEQ_TARGET_DEPENDENCIES})
target_link_libraries(mmseqs_native PRIVATE ${MMSEQ_TARGET_LINKS})

# Pedantic
if (MSVC)
    # warning level 4 and all warnings as errors
    target_compile_options(mmseqs_native PRIVATE /W4 /WX)
else()
    # lots of warnings and all warnings as errors
    target_compile_options(mmseqs_native PRIVATE -Wall -Wextra -pedantic)
endif()

# needed for concat.h
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>

        int main() {
          FILE* tmpf = tmpfile();
          int input_desc = fileno(tmpf);
          int test = posix_fadvise(input_desc, 0, 0, POSIX_FADV_SEQUENTIAL);
          fclose(tmpf);
          return 0;
        }"
        HAVE_POSIX_FADVISE)
if (HAVE_POSIX_FADVISE)
    target_compile_definitions(mmseqs_native PUBLIC -DHAVE_POSIX_FADVISE=1)
endif ()

check_cxx_source_compiles("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>
        #include <sys/mman.h>

        int main() {
          FILE* tmpf = tmpfile();
          void *tmp = mmap(NULL, 32, PROT_READ, MAP_SHARED, fileno(tmpf), 0);
          int test = posix_madvise(tmp, 32, POSIX_MADV_SEQUENTIAL);
          fclose(tmpf);
          return 0;
        }"
        HAVE_POSIX_MADVISE)
if (HAVE_POSIX_MADVISE)
    target_compile_definitions(mmseqs_native PUBLIC -DHAVE_POSIX_MADVISE=1)
endif ()

if (NOT DISABLE_IPS4O)
    find_package(Atomic)
    if (ATOMIC_FOUND)
        target_link_libraries(mmseqs_native PRIVATE ${ATOMIC_LIBRARIES})
        target_compile_definitions(mmseqs_native PUBLIC -DENABLE_IPS4O=1)
        message("-- IPS4O sorting works")
    else ()
        message("-- OMPTL sorting fallback")
    endif ()
else ()
    message("-- OMPTL sorting fallback")
endif ()

target_link_libraries(mmseqs_native PRIVATE tinyexpr libzstd_static microtar)
if (CYGWIN)
    target_link_libraries(mmseqs_native PRIVATE nedmalloc)
endif ()

find_package(ZLIB QUIET)
if (ZLIB_FOUND)
    message("-- Found ZLIB")
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(OLD_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS -lz)
    set(CMAKE_REQUIRED_INCLUDES ${ZLIB_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${ZLIB_LIBRARIES})
    check_cxx_source_compiles("
        #include <zlib.h>
        int main() { gzFile file; return 0; }"
        HAVE_ZLIB_CHECK)
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OLD_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${OLD_CMAKE_REQUIRED_LIBRARIES})
    if(HAVE_ZLIB_CHECK)
        message("-- ZLIB works")
        target_include_directories(mmseqs_native PUBLIC ${ZLIB_INCLUDE_DIRS})
        target_compile_definitions(mmseqs_native PUBLIC -DHAVE_ZLIB=1)
        target_link_libraries(mmseqs_native PRIVATE ${ZLIB_LIBRARIES})
    else ()
        message("-- ZLIB does not work")
    endif()
else ()
    message("-- Could not find ZLIB")
endif ()

find_package(BZip2 QUIET)
if (BZIP2_FOUND)
    message("-- Found BZLIB")
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(OLD_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS -lbz2)
    set(CMAKE_REQUIRED_INCLUDES ${BZIP_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${BZIP2_LIBRARIES})
    check_cxx_source_compiles("
        #include <bzlib.h>
        int main() { bz_stream stream; return 0; }"
        HAVE_BZLIB_CHECK)
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OLD_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${OLD_CMAKE_REQUIRED_LIBRARIES})
    if(HAVE_BZLIB_CHECK)
        message("-- BZLIB works")
        target_include_directories(mmseqs_native PUBLIC ${BZIP_INCLUDE_DIRS})
        target_compile_definitions(mmseqs_native PUBLIC -DHAVE_BZLIB=1)
        target_link_libraries(mmseqs_native PRIVATE ${BZIP2_LIBRARIES})
    else ()
        message("-- BZLIB does not work")
    endif()
else ()
    message("-- Could not find BZLIB")
endif ()

# MPI
if (HAVE_MPI)
    find_package(MPI REQUIRED)
    if (MPI_FOUND)
        message("-- Found MPI")
        target_include_directories(mmseqs_native PUBLIC ${MPI_INCLUDE_PATH})
        #Hack (OMPI_SKIP_MPICXX=1): https://github.com/open-mpi/ompi/issues/5157#issuecomment-388495496
        target_compile_definitions(mmseqs_native PUBLIC -DHAVE_MPI=1 -DOMPI_SKIP_MPICXX=1)
        target_link_libraries(mmseqs_native ${MPI_LIBRARIES})
        append_target_property(mmseqs_native COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
        append_target_property(mmseqs_native LINK_FLAGS ${MPI_LINK_FLAGS})
    endif ()
endif ()

find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    message("-- Found OpenMP")
    target_compile_definitions(mmseqs_native PUBLIC -DOPENMP=1)
    # For GCC we dont want to do this since it breaks macOS static builds
    # It will link libgomp.a internally (through -fopenmp I guess)
    # and also link libgomp.dylib thus breaking static builds
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_link_libraries(mmseqs_native PRIVATE ${OpenMP_CXX_LIBRARIES})
    endif()
    append_target_property(mmseqs_native COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    append_target_property(mmseqs_native LINK_FLAGS ${OpenMP_CXX_FLAGS})
elseif (REQUIRE_OPENMP)
    message("-- Could not find OpenMP.")
endif ()

if (HAVE_GPROF)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-pg GPROF_FOUND)
    if (GPROF_FOUND)
        append_target_property(mmseqs_native COMPILE_FLAGS -pg)
        append_target_property(mmseqs_native LINK_FLAGS -pg)
    else ()
        message(FATAL_ERROR "-- Could not find GPROF")
    endif ()
endif ()

###

# add_executable(state_loader src/state/state_loader.cpp)
# target_include_directories(state_loader PUBLIC ${MMSEQ_TARGET_INCLUDES})
# add_dependencies(state_loader ${MMSEQ_TARGET_DEPENDENCIES})
# target_link_libraries(state_loader PRIVATE ${MMSEQ_TARGET_LINKS})


