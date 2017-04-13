cmake_minimum_required(VERSION 3.6)
project(LDA)

set(CMAKE_CXX_STANDARD 11)

include_directories("/Users/lixupeng/tbb/include")
link_directories("/Users/lixupeng/tbb/lib")

set(SOURCE_FILES main.cpp maps.h types.h FTreeLDA.h TraverseHashMap.h FTree.h BinarySearch.h)
add_executable(LDA ${SOURCE_FILES})
target_link_libraries(LDA tbb tbbmalloc)
