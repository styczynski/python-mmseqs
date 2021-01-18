cmake_minimum_required(VERSION 2.8.12)
project(unafold_python_native)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(PKGDATADIR ${CMAKE_INSTALL_FULL_DATADIR}/pyunafold)
add_definitions(-DPKGDATADIR="${PKGDATADIR}" -D__GNU_LIBRARY__ -DEXIT_SUCCESS=0 -DEXIT_FAILURE=1 -DACKAGE_BUGREPORT=" <nowhere>" -DPACKAGE_NAME="pyunafold" -DPACKAGE_VERSION="1.0.0")

SET(SOURCE_DIR "src/unafold_python_native")
# Tell cmake that headers are in alse in source_dir
include_directories(${SOURCE_DIR})
SET(SOURCES "${SOURCE_DIR}/unafold/energy.cpp" "${SOURCE_DIR}/unafold/xmalloc.cpp" "${SOURCE_DIR}/unafold/hybrid-min.cpp")

# Generate python module
add_subdirectory(lib/pybind11)
pybind11_add_module(unafold_python_native ${SOURCES} "${SOURCE_DIR}/bindings.cpp")