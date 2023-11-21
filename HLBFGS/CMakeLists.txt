#cmake_minimum_required(VERSION 2.8.9)
#project(HLBFGS)
#set(CMAKE_CXX_STANDARD 11)

#aux_source_directory(${PROJECT_SOURCE_DIR} HLBFGS_FILES)
#add_library(${PROJECT_NAME} STATIC ${HLBFGS_FILES})

cmake_minimum_required(VERSION 2.8.9)
project(HLBFGS)
set(CMAKE_CXX_STANDARD 11)

#set(LIBRARY_OUTPUT_PATH "build")
aux_source_directory(${PROJECT_SOURCE_DIR} HLBFGS_src)
include_directories(${PROJECT_SOURCE_DIR})
add_library(HLBFGS SHARED ${HLBFGS_src})