project (FEvoR)

cmake_minimum_required (VERSION 2.8.0)
set(CMAKE_MACOSX_RPATH 0)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR})

option (BUILD_SHARED_LIBS "Build shared libraries" ON)

include_directories(src)

set(test_sources
  test/test_crystal.cc
  test/test_distribution.cc
  )

add_library(fevor_common
  src/fevor_crystal.cc
  src/fevor_distribution.cc
  src/vector_tensor_opperations.cc
)

add_executable(singleStep Example/singleStep.cc)
target_link_libraries(singleStep fevor_common)

add_executable(testUnits ${test_sources} test/units.cc)
target_link_libraries(testUnits fevor_common)

# This one compiles, but will crash unless executed from the top
# directory of the source tree (the path to an input file is
# hard-wired).
add_executable(comparisonStep ${test_sources} test/comparisonStep.cc)
target_link_libraries(comparisonStep fevor_common)
