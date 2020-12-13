cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 11)

add_definitions     ( ${EIGEN3_DEFINITIONS} )
include_directories ( ${EIGEN3_INCLUDE_DIRS} )

project(solver.cpp)

find_package(OsqpEigen REQUIRED)
find_package(Eigen3)
# find_package(kindr)

include_directories(
    SYSTEM ${EIGEN3_INCLUDE_DIR} 
    # ${kindr_INCLUDE_DIRS} 
    include
)


add_compile_options(
    -D_GLIBCXX_USE_CXX11_ABI=0
    -DEIGEN_STACK_ALLOCATION_LIMIT=0
)

#solver.cpp
add_executable(
    solver.cpp 
    src/solver.cpp 
    src/ModelPredictiveControlAPI.cpp 
    include/ModelPredictiveControlAPI.h
    src/SerialPort.cpp
    include/SerialPort.h
)

target_link_libraries(solver.cpp OsqpEigen::OsqpEigen)