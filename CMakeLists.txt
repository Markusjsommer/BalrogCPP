cmake_minimum_required(VERSION 3.16)
project(Balrog)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

file(GLOB globbed
        "*.h"
        "*.cpp"
        )

include_directories(include)
include_directories(include/gzip)

#SET(Torch_DIR /home/markus/libtorch/share/cmake/Torch)
#SET(Torch_DIR libtorch/share/cmake/Torch)

find_package(Torch REQUIRED)

add_executable(Balrog main.cpp ${globbed} include/tqdm.h include/FastaReader.cpp include/FastaReader.h include/GeneFinder.cpp include/GeneFinder.h)

target_link_libraries(Balrog ${TORCH_LIBRARIES})
set_property(TARGET Balrog PROPERTY CXX_STANDARD 14)


include(FindPkgConfig)
find_package(zlib REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})

target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARY})

