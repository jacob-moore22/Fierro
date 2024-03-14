cmake_minimum_required(VERSION 3.1.3)


find_package(Kokkos REQUIRED)
add_definitions(-DHAVE_KOKKOS=1)

include_directories(include)
include_directories(src)

set(SRC_Files 
src/boundary.cpp 
src/energy_sgh.cpp
src/eos.cpp
src/force_sgh.cpp
src/geometry.cpp
src/momentum.cpp
src/properties.cpp
src/sgh_solve.cpp
src/time_integration.cpp
src/user_mat.cpp
src/user_mat_init.cpp
src/outputs.cpp)

set(SGH_Solver_SRC src/sgh_solver.cpp )

message("\n ****** ADDING SGH LIBRARY ******** \n ")

add_library(sgh_solver OBJECT ${SRC_Files} ${SGH_Solver_SRC}) # 

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()

target_link_libraries(sgh_solver Kokkos::kokkos)


