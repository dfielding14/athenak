cmake_minimum_required(VERSION 3.16)
project(gotham_pdf_rebuild LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(ATHENAK_SOURCE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." CACHE PATH
    "AthenaK source checkout providing the bundled Kokkos")

find_package(MPI REQUIRED COMPONENTS CXX)

if(Kokkos_ENABLE_CUDA)
  set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "" FORCE)
endif()

add_subdirectory("${ATHENAK_SOURCE_DIR}/kokkos" "${CMAKE_BINARY_DIR}/kokkos"
                 EXCLUDE_FROM_ALL)

add_executable(gotham_pdf_rebuild gotham_pdf_rebuild.cpp)
target_link_libraries(gotham_pdf_rebuild PRIVATE Kokkos::kokkos MPI::MPI_CXX)
