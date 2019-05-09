cmake_minimum_required(VERSION 3.11)
project(NERE_2D)

set(CMAKE_CXX_STANDARD 11)

find_package(Armadillo REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)

include_directories(${ARMADILLO_PATH})
include_directories(${LAPACK_PATH})
include_directories(${EIGEN3_PATH})

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Winline")

if (CMAKE_BUILD_TYPE MATCHES Debug)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O1")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
endif()

add_executable(NERE_2D main.cpp Richards.cpp Richards.h global.h KFun.cpp KFun.h hysteresis.cpp hysteresis.h melt.h melt.cpp)


target_link_libraries(NERE_2D ${ARMADILLO_LIBRARIES} ${LAPACK_LIBRARIES} Eigen3::Eigen)

