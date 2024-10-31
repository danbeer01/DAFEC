program Main
!
! Purpose:
! To read from an input file and then calculate flux values and k_eff for a multi-group, 
! multi-material, fission source problem using the Finite Element Method.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 07/03/2024    Daniel Beer     Original code
!
use C_timer
use m_VTK_Reader
use m_Read_Properties
use m_Solver
use m_Results
use m_Output
use m_VTK_Writer

implicit none

! Variable declarations
type(PropertiesType) :: Properties
type(NType)          :: N
type(ResultsType)    :: Results
type(t_timer)        :: timer
type(Mesh)           :: this

! Start timer
call timer%start_timer()
call timer%Startdate()

! Open VTK file
open(unit=1, file='vtk/O1_line_3.vtk', status='old', action='read')
! open(unit=1, file='vtk/reed_2.vtk', status='old', action='read')
! open(unit=1, file='vtk/3_Region_Sphere.vtk', status='old', action='read')

! Read mesh from VTK file
call this%read_VTK_file()

! Open input file for material and other data data
open (unit=1, file='inputs/homo.txt', status='old', action='read')
! open (unit=1, file='inputs/reed.txt', status='old', action='read')
! open (unit=1, file='inputs/sood_3.txt', status='old', action='read')
! open (unit=1, file='inputs/3_Region_Sphere.txt', status='old', action='read')

! Read properties from input file
call Read_Properties(Properties, N, this)

! Solve for flux and k_eff
call Solver(Properties, Results, N, this)

! Output flux and k_eff
! call neutron_balance(Properties, Results, N, this)
! call reaction_rate(Properties, Results, N, this)
call print_flux(Properties, Results, N, this)
if(Properties%Case == 1) call Results%print_keff()

! End timer
call timer%elapsed_time()
call timer%elapseddate()

! Output VTK file
call Output_VTK(this, Results, N)

end program Main