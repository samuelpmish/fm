enable_testing()

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

file(GLOB cpp_tests ${PROJECT_SOURCE_DIR}/tests/*.cpp)

#add_library(testlib ${PROJECT_SOURCE_DIR}/tests/pch.cpp)
#target_precompile_headers(testlib PRIVATE ${PROJECT_SOURCE_DIR}/tests/pch.cpp)
#target_link_libraries(testlib 

foreach(filename ${cpp_tests})		

  get_filename_component(testname ${filename} NAME_WE)

  add_executable(${testname} ${filename})
  target_link_libraries(${testname} PUBLIC fm GTest::gtest_main)

  add_test(${testname} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${testname})		

endforeach(filename ${cpp_tests})

if (FM_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()

  enable_language(CUDA)

  file(GLOB cu_tests ${PROJECT_SOURCE_DIR}/tests/*.cu)

  foreach(filename ${cu_tests})		

    get_filename_component(testname ${filename} NAME_WE)

    add_executable(cuda_${testname} ${filename})
    target_link_libraries(cuda_${testname} PUBLIC fm GTest::gtest_main)

    add_test(cuda_${testname} ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cuda_${testname})		

  endforeach(filename ${cu_tests})

endif()