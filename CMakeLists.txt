cmake_minimum_required(VERSION 2.8.9)
project(balance_solver)
set(CMAKE_CXX_STANDARD 11)



# 设置头文件&相关库
include_directories(include)
include_directories(eigen3) 
find_package(OpenMP REQUIRED)
#find_package(Torch REQUIRED)

#add_subdirectory(HLBFGS)
add_subdirectory(pybind11)
#include_directories(HLBFGS)
# 模块
pybind11_add_module(balance_solver ./src/balance_solver.cpp)
#target_link_libraries(${PROJECT_NAME} HLBFGS)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)