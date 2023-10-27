cmake_minimum_required(VERSION 2.8.9)
project(pybind_eigen_test)
set(CMAKE_CXX_STANDARD 11)

# 设置头文件&相关库
include_directories(include)
include_directories(eigen3) 

add_subdirectory(pybind11)

# 模块
pybind11_add_module(pd_cpp ./src/pd.cpp)