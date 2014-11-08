cmake_minimum_required (VERSION 2.8.0)

project (FEvoR)
set (${PROJECT_NAME}_VERSION_MAJOR 0)
set (${PROJECT_NAME}_VERSION_MINOR 1)

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
  src/Faddeeva.cc
)

add_executable(singleStep Example/singleStep.cc)
target_link_libraries(singleStep FEvoR)

add_executable(enhancementFactor Example/enhancementFactor.cc)
target_link_libraries(enhancementFactor FEvoR)


add_executable(testUnits ${test_sources} test/units.cc)
target_link_libraries(testUnits FEvoR)

add_executable(comparisonStep test/comparisonStep.cc)
target_link_libraries(comparisonStep FEvoR)

# build package config file
set(INCLUDE_INSTALL_DIR include CACHE PATH "headers directory" )
set(LIB_INSTALL_DIR lib CACHE PATH "library directory" )

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${PROJECT_NAME}Config.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${LIB_INSTALL_DIR}/${PROJECT_NAME}
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
    COMPATIBILITY AnyNewerVersion )

# install the FEvoR library
install(TARGETS FEvoR
        LIBRARY DESTINATION lib/${PROJECT_NAME})

install(FILES 
  src/fevor_distribution.hh
  src/fevor_crystal.hh
  src/vector_tensor_operations.hh
  src/Faddeeva.hh
  DESTINATION include/${PROJECT_NAME})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake 
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION lib/${PROJECT_NAME})

