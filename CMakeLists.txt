cmake_minimum_required(VERSION 3.11)
project(Dyn_2D)
set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
find_package(lapack REQUIRED)
#find_package(Eigen3 REQUIRED)

include_directories(${ARMADILLO_PATH})
include_directories(${LAPACK_PATH})
#include_directories(${EIGEN3_PATH})

set(CMAKE_CXX_COMPILER /usr/local/bin/g++-8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -Winline -O2")

add_executable(NERE_2D main.cpp Richards.cpp Richards.h global.h KFun.cpp KFun.h hysteresis.cpp hysteresis.h)


target_link_libraries(NERE_2D ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES})

