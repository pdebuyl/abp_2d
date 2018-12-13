cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT "Release")

Project(abp_2d Fortran)
find_package(OpenMP)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

add_library(abp_model src/common.f90 src/abp_model_m.f90)

add_executable(test_min_dist test/test_min_dist.f90)
target_link_libraries(test_min_dist abp_model)

add_executable(test_noise_term test/test_noise_term.f90)
target_link_libraries(test_noise_term abp_model)

add_executable(run_abp_probe src/run_abp_probe.f90)
target_link_libraries(run_abp_probe abp_model)

