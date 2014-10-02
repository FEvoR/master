cmake_minimum_required (VERSION 2.8.0)

project (FEvoR)

# TODO: need to set version numbers. 
# TODO: Build and Export the FEvoRConfig.cmake file automatically
#SET(${PROJECT_NAME}_VERSION_MAJOR 0)
#SET(${PROJECT_NAME}_VERSION_MINOR 1)

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

add_library(FEvoR
  src/fevor_crystal.cc
  src/fevor_distribution.cc
  src/vector_tensor_operations.cc
)

add_executable(singleStep Example/singleStep.cc)
target_link_libraries(singleStep FEvoR)

add_executable(testUnits ${test_sources} test/units.cc)
target_link_libraries(testUnits FEvoR)

# This one compiles, but will crash unless executed from the top
# directory of the source tree (the path to an input file is
# hard-wired).
add_executable(comparisonStep ${test_sources} test/comparisonStep.cc)
target_link_libraries(comparisonStep FEvoR)

# install the FEvoR library
install(TARGETS FEvoR
        LIBRARY DESTINATION lib/${CMAKE_PROJECT_NAME})

install(FILES 
  src/fevor_distribution.hh
  src/fevor_crystal.hh
  src/vector_tensor_operations.hh
  DESTINATION include/${CMAKE_PROJECT_NAME})

install(FILES FEvoRConfig.cmake DESTINATION lib/${CMAKE_PROJECT_NAME})

