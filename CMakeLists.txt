cmake_minimum_required(VERSION 3.16)
project(balrog)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

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

if (APPLE)
  set(CMAKE_INSTALL_RPATH "@executable_path/../lib/python3.7/site-packages/torch/lib")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/python3.7/site-packages/torch/lib")
endif()

find_package(Torch REQUIRED)

add_executable(balrog main.cpp ${globbed} include/tqdm.h include/FastaReader.cpp include/FastaReader.h include/GeneFinder.cpp include/GeneFinder.h)

target_link_libraries(balrog ${TORCH_LIBRARIES})
set_property(TARGET balrog PROPERTY CXX_STANDARD 14)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( balrog ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )


include(CMakeRC)
cmrc_add_resource_library(
        cmakeresources
        NAMESPACE cmakeresources
        WHENCE data
        data/gene_model_v1.0.pt
        data/TIS_model_v1.0.pt
        data/reference_genes.tar.gz
        )

target_link_libraries(${PROJECT_NAME} cmakeresources)

install(
  TARGETS balrog
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
)
