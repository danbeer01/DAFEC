cmake_minimum_required(VERSION 3.6) 

# Set our project name and Fortran as only language we compile 
# (C and C++ are also available) and should be added if we interop 
project(neutron_code VERSION 0.1 
    DESCRIPTION "Daniels"
    LANGUAGES Fortran C
    )

# we default to Release build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
#     set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall")
#     set(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback")
#     set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHOST")
# endif()

if(APPLE)
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -ffree-line-length-none ")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-ffree-line-length-800 -O0 -fcheck=bounds -DDEBUG -pg -fopenmp ")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native -DNDEBUG -fopenmp ")
    endif()
else()
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
        set(CMAKE_Fortran_FLAGS         "${CMAKE_Fortran_FLAGS} -Wall -ffree-line-length-none ")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-ffree-line-length-800 -O0 -g -pg -fcheck=bounds -fdec -DDEBUG -pg -fopenmp ")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -march=native -DNDEBUG -fopenmp ")
    endif()
endif()


set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/..)

# set(CMAKE_OSX_ARCHITECTURES x86_64)

# set(CMAKE_VERBOSE_MAKEFILE ON)

# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcrt1.o")

# EXTERNAL PACKAGES
# Fortran Standard Library:

find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)

#find_package(fortran_stdlib REQUIRED)

#FIND_PACKAGE(PETSc REQUIRED)

# To use this please have stdlib built with CMake already 

#NB this should be the actual build e.g. the ARCH folder I think
# (I am using a ./configure --prefix <....> build)
set(MY_PETSC_PATH "/home/daniel/Documents/apps") 
set(CMAKE_PREFIX_PATH ${MY_PETSC_PATH}) 

# the `pkg_check_modules` function is created with this call
find_package(PkgConfig REQUIRED)
pkg_check_modules(PETSC REQUIRED PETSc)

# Add include directories
INCLUDE_DIRECTORIES(${PETSC_INCLUDE_DIRS})
message(${PETSC_INCLUDE_DIRS})

# Add our source files
# All .f90 files in lib/ will be added to the library as sources
# Below command finds all .f90 and adds them to the ${SOURCES} variable

# file(GLOB SOURCES "stdlib/*.f90" "stdlib/*.F90" main.f90)
# add_executable(output ${SOURCES})

add_library(src
    src/C_timer.f90
    src/m_Constants_Mod.f90
    src/m_Determine_Solve_Method.f90
    src/m_Gauss_Points.f90
    src/m_VTK_Reader.f90

    src/Transport/m_Read_Properties.f90
    src/Transport/m_Results.f90
    src/Transport/m_Small_Matrix_Solver.f90
    src/Transport/m_Normalisation.f90
    src/Transport/m_Create_Hexahedral_Shape_Functions.f90
    src/Transport/m_Create_Tetrahedral_Shape_Functions.f90
    src/Transport/m_Create_Prismatic_Shape_Functions.f90
    src/Transport/m_Create_Pyramidal_Shape_Functions.f90
    src/Transport/m_Create_Quadrilateral_Shape_Functions.f90
    src/Transport/m_Create_Triangular_Shape_Functions.f90
    src/Transport/m_Create_Shape_Functions.f90
    src/Transport/m_Construct_Matrix_1D.f90
    src/Transport/m_Construct_Matrix_2D.f90
    src/Transport/m_Construct_Matrix_3D.f90
    src/Transport/m_Calculate_mu_w.f90
    src/Transport/m_Calculate_Curvilinear_Normals.f90
    src/Transport/m_Sweep_Order_1D.f90
    src/Transport/m_Sweep_Order_2D.f90
    src/Transport/m_Sweep_Order_3D.f90
    src/Transport/m_Boundary_Conditions.f90
    src/Transport/m_Cylindrical.f90
    src/Transport/m_Spherical.f90
    src/Transport/m_RZ.f90
    src/Transport/m_Solver.f90
    src/Transport/m_Output.f90
    src/Transport/m_VTK_Writer.f90

    src/Diffusion/m_Read_Properties_D.f90
    src/Diffusion/m_Boundary_Conditions_D.f90
    src/Diffusion/m_Results_D.f90
    src/Diffusion/m_PETSc.F90
    src/Diffusion/m_Normalisation_D.f90
    src/Diffusion/m_Create_Hexahedral_Shape_Functions_D.f90
    src/Diffusion/m_Create_Tetrahedral_Shape_Functions_D.f90
    src/Diffusion/m_Create_Prismatic_Shape_Functions_D.f90
    src/Diffusion/m_Create_Pyramidal_Shape_Functions_D.f90
    src/Diffusion/m_Create_Quadrilateral_Shape_Functions_D.f90
    src/Diffusion/m_Create_Triangular_Shape_Functions_D.f90
    src/Diffusion/m_Create_Shape_Functions_D.f90
    src/Diffusion/m_Construct_Matrix_D.f90
    src/Diffusion/m_Solver_D.f90
    src/Diffusion/m_Output_D.f90
    src/Diffusion/m_VTK_Writer_D.f90
    
    )

# add_subdirectory(stdlib)

add_executable(out main.f90)

target_link_libraries(out src ${PETSC_LINK_LIBRARIES})

# target_link_libraries(output LAPACK::LAPACK BLAS::BLAS)

#TARGET_LINK_LIBRARIES(out ${PETSC_LINK_LIBRARIES} fortran_stdlib::fortran_stdlib)
