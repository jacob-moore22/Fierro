project (FierroSGH)


cmake_minimum_required(VERSION 3.1.3)

# include_directories(../../../../single-node/install/kokkos-serial/kokkos/include)
# include_directories(../../../../single-node/install/kokkos-serial/kokkos/bin)

# include_directories(../../../../single-node/install/kokkos-serial/kokkos/lib64/cmake/Kokkos)
# add_subdirectory(../../../../single-node/install/kokkos-serial/kokkos/lib64/cmake/Kokkos)


include_directories(../../../../lib/Elements/matar/src)
include_directories(../../../../lib/Elements/elements)
include_directories(../../../../single-node/install/kokkos-serial/kokkos/include)
include_directories(../../common)


# add_subdirectory(../../../../lib/Elements cbin)

include_directories(../../../../lib/Elements)
add_subdirectory(../../../../lib/Elements cbin)
include_directories(../../../../lib/Elements/common)

include_directories(../../)

# include_directories(../../../../lib/Elements/matar)
# add_subdirectory(../../../../lib/Elements/matar cbin)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /var/tmp/repos/Fierro/fork/Fierro/single-node/install/kokkos-serial/kokkos/lib64/cmake/Kokkos/KokkosConfig.cmake)


set(Kokkos_DIR /var/tmp/repos/Fierro/fork/Fierro/single-node/install/kokkos-serial/kokkos/lib64/cmake/Kokkos)
include(/var/tmp/repos/Fierro/fork/Fierro/single-node/install/kokkos-serial/kokkos/lib64/cmake/Kokkos/KokkosConfig.cmake)


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
src/user_mat_init.cpp)

set(SGH_Module_SRC src/sgh_solver.cpp )


# add_executable(FierroSGH driver.cpp read_mesh.cpp geometry.cpp boundary.cpp input.cpp eos.cpp setup.cpp outputs.cpp sgh_solve.cpp time_integration.cpp momentum.cpp force_sgh.cpp 
# energy_sgh.cpp properties.cpp user_mat.cpp user_mat_init.cpp)

if (CUDA)
  add_definitions(-DHAVE_CUDA=1)
elseif (HIP)
  add_definitions(-DHAVE_HIP=1)
elseif (OPENMP)
  add_definitions(-DHAVE_OPENMP=1)
elseif (THREADS)
  add_definitions(-DHAVE_THREADS=1)
endif()



add_library(sgh_solver OBJECT ${SRC_Files} ${SGH_Module_SRC})


# target_link_libraries(FierroSGH matar Kokkos::kokkos)

