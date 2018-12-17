cmake_minimum_required(VERSION 3.1)
set(CMAKE_BUILD_TYPE_INIT "Release")

Project(abp_2d C Fortran)
find_package(OpenMP)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

include_directories("include")

add_library(abp_model_common src/common.f90)

add_library(abp_model src/config_m.f90  src/abp_model_m.f90)
target_link_libraries(abp_model abp_model_common)

add_library(threefry_m src/threefry_c_module.c src/threefry_m.f90)
set_property(TARGET threefry_m PROPERTY C_STANDARD 99)
target_link_libraries(threefry_m abp_model_common)

add_executable(test_min_dist test/test_min_dist.f90)
target_link_libraries(test_min_dist abp_model)

add_executable(test_noise_term test/test_noise_term.f90)
target_link_libraries(test_noise_term abp_model)

add_executable(run_abp_probe src/run_abp_probe.f90)
target_link_libraries(run_abp_probe abp_model)

add_executable(test_threefry_0 test/test_threefry_0.f90)
target_link_libraries(test_threefry_0 threefry_m)

add_executable(test_threefry_1 test/test_threefry_1.f90)
target_link_libraries(test_threefry_1 threefry_m)
