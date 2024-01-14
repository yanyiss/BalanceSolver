cmake_minimum_required(VERSION 3.11)
project(balance_solver)
set(CMAKE_CXX_STANDARD 11)



# 设置头文件&相关库
include_directories(include)
include_directories(eigen3) 
#find_package(OpenMP REQUIRED)
#set(Torch_DIR "/home/yanyisheshou/anaconda3/envs/code/lib/python3.9/site-packages/torch/share/cmake/Torch")
#find_package(Torch REQUIRED)

#add_subdirectory(HLBFGS)
add_subdirectory(pybind11)
#include_directories(HLBFGS)
# 模块
pybind11_add_module(balance_solver ./src/balance_solver.cpp)
#target_link_libraries(${PROJECT_NAME} "${TORCH_LIBRARIES}")
#target_link_libraries(${PROJECT_NAME} HLBFGS)
#target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)