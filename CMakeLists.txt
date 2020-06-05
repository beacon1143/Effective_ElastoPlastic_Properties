cmake_minimum_required(VERSION 3.1)

project(EffProps)

set(EFF_PROPS_SRC main.cpp SigmaCalc.cpp InputData.cpp)
set(EFF_PROPS_HDR SigmaCalc.hpp InputData.hpp)

add_executable(${PROJECT_NAME} ${EFF_PROPS_SRC} ${EFF_PROPS_HDR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")