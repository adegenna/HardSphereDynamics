cmake_minimum_required(VERSION 2.6)
project(lagrangiansolver)

find_package(Boost COMPONENTS program_options REQUIRED)
find_package(GTest REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Threads)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_definitions( -DSRCDIR=\"${CMAKE_SOURCE_DIR}/\")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Werror=reorder")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -DBOOST_LOG_DYN_LINK -DNDEBUG")

include_directories( ${CMAKE_SOURCE_DIR}
                     ${Boost_INCLUDE_DIRS}
                     ${GTEST_INCLUDE_DIRS}
                     ${EIGEN3_INCLUDE_DIR}
)

add_executable( testall
  tests/tests.cpp
  src/Inputfile.hpp
  src/LagrangianState.cpp
  src/ParticlePhysics.cpp
  src/TimeIntegration.cpp
  tests/inputs_test.cpp
  tests/lagrangian_test.cpp
  tests/physics_test.cpp
)

add_executable( solver
  src/driver.cpp
  src/Inputfile.hpp
  src/LagrangianState.cpp
  src/ParticlePhysics.cpp
  src/TimeIntegration.cpp
  src/options.hpp
  src/options_parser.hpp
)


target_link_libraries( solver
  ${Boost_LIBRARIES}   
)

target_link_libraries(testall
  ${GTEST_LIBRARIES}
  ${Boost_LIBRARIES}
  pthread
)


enable_testing()

set_target_properties( testall PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

add_test(
  NAME
  testall
  COMMAND
  testall
)
