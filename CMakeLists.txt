cmake_minimum_required(VERSION 3.20)
project(ahc_cpp)

set(CMAKE_CXX_STANDARD 17)

set(AHC_SRC include/ahc/ahc.h src/agglomerative_clustering.h src/agglomerative_clustering.cpp src/util.h src/ahc_impl.cpp)

add_library(ahc STATIC ${AHC_SRC})
target_include_directories(ahc PUBLIC include)