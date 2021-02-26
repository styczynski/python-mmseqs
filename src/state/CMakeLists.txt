cmake_minimum_required(VERSION 2.8.12)
project(state_loader)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake")

set(CMAKE_BUILD_TYPE Debug)
if (MSVC)
    add_compile_options(/DEBUG:FULL)
else ()
    add_compile_options(-g)
endif ()

# cppfs
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/cppfs build_cppfs)

# spdlog
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../../lib/spdlog build_spdlog)

add_executable(state_loader state_loader.c)

target_include_directories(state_loader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/toml11)
target_include_directories(state_loader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/spdlog/include)
#target_include_directories(state_loader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../lib/cppfs/include)
target_include_directories(state_loader PUBLIC .)

target_link_libraries(state_loader PRIVATE spdlog)
#target_link_libraries(state_loader PRIVATE cppfs)
