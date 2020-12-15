cmake_minimum_required(VERSION 3.16)
project(balrog)

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

execute_process(
  COMMAND python -c "import torch; import os; print(os.path.dirname(torch.__file__), end='')"
  OUTPUT_VARIABLE TorchPath
)
list(APPEND CMAKE_PREFIX_PATH ${TorchPath})

set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#set(CMAKE_INSTALL_RPATH "$ORIGIN/../torch/lib")

set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/python3.7/site-packages/torch/lib")

#execute_process(
#  COMMAND python -c "import torch; import os; print(os.path.join(os.path.dirname(torch.__file__), 'lib'), end='')"
#  OUTPUT_VARIABLE LD_LIBRARY_PATH
#)

find_package(Torch REQUIRED)

add_executable(balrog main.cpp ${globbed} include/tqdm.h include/FastaReader.cpp include/FastaReader.h include/GeneFinder.cpp include/GeneFinder.h)

target_link_libraries(balrog ${TORCH_LIBRARIES})
set_property(TARGET balrog PROPERTY CXX_STANDARD 14)



find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( balrog ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )


#target_link_libraries(${PROJECT_NAME} -lz)


#include(FindPkgConfig)
#find_package(zlib REQUIRED)
#include_directories(${ZLIB_INCLUDE_DIR})

#target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARY})

install(
  TARGETS balrog
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
