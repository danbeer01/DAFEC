program Main
!
! Purpose:
! To read from an input file and then calculate flux values and k_eff for a multi-group, 
! multi-material, fission source problem using the Finite Element Method.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
use C_timer
use m_VTK_Reader
use Determine_Solve_Method
use m_Read_Properties
use m_Read_Properties_D
use m_Solver
use m_Solver_D
use m_Results
use m_Results_D
use m_Output
use m_Output_D
use m_VTK_Writer
use m_VTK_Writer_D

implicit none

! Variable declarations
type(PropertiesType)  :: Properties
type(NType)           :: N
type(PropertiesTypeD) :: PropertiesD
type(NTypeD)          :: N_D
type(ResultsType)     :: Results
type(ResultsTypeD)    :: ResultsD
type(t_timer)         :: timer
type(Mesh)            :: vtk_mesh

logical :: Solve_Transport_Method

! Start timer
call timer%start_timer()
call timer%Startdate()

! Open VTK file
! open(unit=1, file='vtk/1D/O1_line_3.vtk', status='old', action='read')
! open(unit=1, file='vtk/1D/reed_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/1D/3_Region_Sphere.vtk', status='old', action='read')

! open(unit=1, file='vtk/2D/Squares_3.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/Tri_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/quad.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O1_tri_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/tri_linear.vtk', status='old', action='read')
open(unit=1, file='vtk/2D/quad_linear.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O2_tri.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O3_tri.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O4_quad.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O5_quad_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/O8_quad.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/C5G2_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/Ackroyd.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/Cylinder.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/LWR.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/mixed_mesh.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/C5G7_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/C5G7vol_1_2_f.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/Curvilinear.vtk', status='old', action='read')
! open(unit=1, file='vtk/2D/pincell_1_2.vtk', status='old', action='read')

! open(unit=1, file='vtk/3D/Cube_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O1_cube.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O1_tet_5.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O2_tet.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O2_hex.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O2_cube.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O3_cube_3.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/Pyr_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O2_Pyr.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/pyramid_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/Pris.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/O2_Pris.vtk', status='old', action='read')
! open(unit=1, file='vtk/3D/Takeda_2.vtk', status='old', action='read')

! Read mesh from VTK file
call vtk_mesh%read_VTK_file()

! Open input file for material and other data data
open(unit=1, file='inputs/homo.txt', status='old', action='read')
! open(unit=1, file='inputs/hetero.txt', status='old', action='read')

! open(unit=1, file='inputs/Benchmarks/reed.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/sood_3.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/3_Region_Sphere.txt', status='old', action='read')

! open(unit=1, file='inputs/Benchmarks/c5g2.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/Ackroyd.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/cylinder.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/LWR.txt', status='old', action='read')
! open(unit=1, file='inputs/Benchmarks/c5g7_3.txt', status='old', action='read')

! open(unit=1, file='inputs/Benchmarks/Takeda.txt', status='old', action='read')

! Determine if problem is transport or diffusion
call Determine_Solve(Solve_Transport_Method)

! Read properties from input file
if(Solve_Transport_Method) call Read_Properties(Properties, N, vtk_mesh)
if(.not. Solve_Transport_Method) call Read_Properties_Diffusion(PropertiesD, N_D, vtk_mesh)

! Solve for flux and k_eff
if(Solve_Transport_Method) call Solver(Properties, Results, N)
if(.not. Solve_Transport_Method) call Solver_Diffusion(PropertiesD, ResultsD, N_D, vtk_mesh)

! Output flux and k_eff
if(Solve_Transport_Method) then
    call neutron_balance(Properties, Results, N)
    call reaction_rate(Properties, Results, N, vtk_mesh)
    call print_flux(Properties, Results, N, vtk_mesh)
else
    call neutron_balance_diffusion(PropertiesD, ResultsD, N_D, vtk_mesh)
    call reaction_rate_diffusion(PropertiesD, ResultsD, N_D, vtk_mesh)
    call print_flux_diffusion(PropertiesD, ResultsD, N_D, vtk_mesh)
end if

if(Properties%Case == 1) call Results%print_keff()
if(PropertiesD%Case == 1) call ResultsD%print_keff_diffusion()

! End timer
call timer%elapsed_time()
call timer%elapseddate()

! Output VTK file
if(Solve_Transport_Method) then
    ! if(N%Degree <= 2) call Output_Discontinuous_VTK(Properties, vtk_mesh, N)
    if(N%Degree <= 2) call Output_VTK(vtk_mesh, Results, N)
    if(N%Degree > 2 .and. N%D == 2) call Output_Higher_Order_VTK(Properties, vtk_mesh, Results, N)
else
    if(N_D%Degree > 2 .and. N_D%D == 2) then
        call Output_Higher_Order_VTK_D(PropertiesD, vtk_mesh, ResultsD, N_D)
    else
        call Output_VTK_D(vtk_mesh, ResultsD, N_D)
    end if
end if

end program Main