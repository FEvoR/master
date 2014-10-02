cmake_minimum_required (VERSION 2.8.0)

project (FEvoR)

set(CMAKE_MACOSX_RPATH 0)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# FEvoR requires the c++11 flag
add_definitions(-std=c++11)

option (BUILD_SHARED_LIBS "Build shared libraries" ON)

include_directories(src)

set(test_sources
  test/test_crystal.cc
  test/test_distribution.cc
  )

add_library(fevor
  src/fevor_crystal.cc
  src/fevor_distribution.cc
  src/vector_tensor_operations.cc
)

add_executable(singleStep Example/singleStep.cc)
target_link_libraries(singleStep fevor)

add_executable(testUnits ${test_sources} test/units.cc)
target_link_libraries(testUnits fevor)

# This one compiles, but will crash unless executed from the top
# directory of the source tree (the path to an input file is
# hard-wired).
add_executable(comparisonStep ${test_sources} test/comparisonStep.cc)
target_link_libraries(comparisonStep fevor)

# install the FEvoR library
install(TARGETS fevor
        LIBRARY DESTINATION lib/${CMAKE_PROJECT_NAME})

install(FILES 
  src/fevor_distribution.hh
  src/fevor_crystal.hh
  src/vector_tensor_operations.hh
  DESTINATION include/${CMAKE_PROJECT_NAME})
