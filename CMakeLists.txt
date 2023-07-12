cmake_minimum_required(VERSION 3.16)

project(fm_proj LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

file(GLOB fm_headers ${PROJECT_SOURCE_DIR}/include/*.hpp)
add_library(fm INTERFACE ${fm_headers})
target_include_directories(fm INTERFACE ${PROJECT_SOURCE_DIR}/include)

if (FM_ENABLE_TESTING)
  add_subdirectory(tests)
endif()
