
cmake_minimum_required(VERSION 3.1)

project (GF2 Fortran C)

enable_language (Fortran C)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#include(FortranCInterface)

set(HDF5_USE_STATIC_LIBRARIES "TRUE")

find_package(HDF5 COMPONENTS C HL Fortran Fortran_HL REQUIRED)

#include(FortranCInterface)
find_package(MPI)
find_package(BLAS)
find_package(LAPACK)


include_directories(${MPI_Fortran_INCLUDE_PATH})
include_directories(${HDF5_INCLUDE_DIRS})
include_directories("${HDF5_INCLUDE_DIRS}/static")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -g -DNDEBUG") # -mkl -assume nounderscore")

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -ffree-line-length-none -fno-underscoring") # -mkl -assume nounderscore")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -mkl") # -assume nounderscore")
endif()

set(SOURCES ${GF2_SOURCE_DIR}/main.f90)

add_subdirectory(${GF2_SOURCE_DIR}/src)
include_directories(${GF2_SOURCE_DIR}/src)
include_directories(${GF2_BINARY_DIR}/src)

add_executable(GF2 ${SOURCES})

target_link_libraries(GF2 ${LAPACK_LIBRARIES}
                          ${BLAS_LIBRARIES}
                          ${MPI_Fortran_LIBRARIES}
                          ${HDF5_C_LIBRARIES}
                          ${HDF5_Fortran_LIBRARIES}
                          ${HDF5_HL_LIBRARIES}
                          ${HDF5_F90_C_LIBRARIES}
                          ${HDF5_HL_F90_C_LIBRARIES}
                          ${HDF5_Fortran_HL_LIBRARIES}
                          gf2-common-lib )
