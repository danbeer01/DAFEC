module m_Read_Properties_D
!
! Purpose:
! To generate the materials for the problem and populate them with data from the input file
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_VTK_Reader

implicit none

type :: NTypeD
    integer :: N
    integer :: D
    integer :: Group
    integer :: Degree
    integer :: Element
    integer :: Material
end type NTypeD

type :: MaterialTypeD
    real(kind=8), dimension(:), allocatable   :: Sigma_a
    real(kind=8), dimension(:), allocatable   :: Sigma_f
    real(kind=8), dimension(:,:), allocatable :: Sigma_s
    real(kind=8), dimension(:), allocatable   :: Sigma_t
    real(kind=8), dimension(:), allocatable   :: Sigma_r
end type MaterialTypeD

type :: ElementTypeD
    integer :: Material
    integer :: Cell_Type
    integer :: Number_of_Nodes
    real(kind=8) :: Volume
    integer, dimension(:), allocatable          :: Cell_Pointers
    real(kind=8), dimension(:,:), allocatable   :: Coordinates
    real(kind=8), dimension(:), allocatable     :: Sigma_a
    real(kind=8), dimension(:), allocatable     :: Sigma_f
    real(kind=8), dimension(:,:), allocatable   :: Sigma_s
    real(kind=8), dimension(:), allocatable     :: Sigma_t
    real(kind=8), dimension(:), allocatable     :: Sigma_r
    real(kind=8), dimension(:,:,:), allocatable :: Jacobian
    real(kind=8), dimension(:,:,:), allocatable :: Inverse_Jacobian
    real(kind=8), dimension(:), allocatable     :: Det_Jacobian
    real(kind=8), dimension(:,:,:), allocatable :: K_Matrix
    real(kind=8), dimension(:,:), allocatable   :: A_Matrix
    real(kind=8), dimension(:,:), allocatable   :: B_Matrix
    real(kind=8), dimension(:,:), allocatable   :: D_Matrix
    real(kind=8), dimension(:,:), allocatable   :: Flux
    real(kind=8), dimension(:,:), allocatable   :: Source
    real(kind=8), dimension(:), allocatable     :: Source_Vector
end type ElementTypeD

type :: PropertiesTypeD
    integer :: LBC,RBC,TBC,BBC,FBC,BaBC
    integer :: Adjoint
    integer :: Case
    integer :: g
    real(kind = 8) :: alpha, Q_s
    real(kind = 8), dimension(:), allocatable     :: Length
    real(kind = 8), dimension(:), allocatable     :: Chi
    real(kind = 8), dimension(:,:), allocatable   :: Isoparametric_Coordinates
    type(ElementTypeD), dimension(:), allocatable  :: Elements
    type(MaterialTypeD), dimension(:), allocatable :: Materials
end type PropertiesTypeD

contains

subroutine Read_Properties_Diffusion(Properties, N, this)
    implicit none

    class(Mesh)                           :: this
    class(PropertiesTypeD), intent(inout) :: Properties
    class(NTypeD), intent(inout)          :: N
    
    integer :: i, j
    real(kind=8), dimension(:), allocatable :: sigma_s_vector

    read(1,*)
    read(1,*) Properties%Case
    read(1,*) N%Material
    read(1,*) N%Group
    read(1,*)
    read(1,*) 
    read(1,*) Properties%g
    read(1,*) Properties%Adjoint
    read(1,*) 
    read(1,*) Properties%LBC
    read(1,*) Properties%RBC
    read(1,*) Properties%TBC
    read(1,*) Properties%BBC
    read(1,*) Properties%FBC
    read(1,*) Properties%BaBC
    read(1,*) Properties%alpha
    read(1,*) Properties%Q_s

    N%Degree = this%N_degree
    
    allocate(Properties%Chi(N%Group))
    
    read(1,*) ! Skip a line to align with input file
    read(1,*) Properties%Chi

    if ((SUM(Properties%Chi) - 1.0_8) > 1e-4) then
        print *, Properties%Chi
        write(*,*) 'Chi values must sum to 1'
        stop
    end if

    N%Element = this%N_cells
    N%D = this%Dimension

    allocate(Properties%Length(N%D))

    do i = 1, N%D
        Properties%Length(i) = MAXVAL(this%Nodes(:,i)) - MINVAL(this%Nodes(:,i))
    end do

    allocate(Properties%Elements(N%Element))
    allocate(Properties%Materials(N%Material))
    allocate(sigma_s_vector(N%Group))

    do i = 1, N%Material
        allocate(Properties%Materials(i)%Sigma_a(N%Group))
        allocate(Properties%Materials(i)%Sigma_f(N%Group))
        allocate(Properties%Materials(i)%Sigma_s(N%Group,N%Group))
        allocate(Properties%Materials(i)%Sigma_T(N%Group))
        allocate(Properties%Materials(i)%Sigma_r(N%Group))
    end do

    do i = 1, N%Element
        allocate(Properties%Elements(i)%Sigma_a(N%Group))
        allocate(Properties%Elements(i)%Sigma_f(N%Group))
        allocate(Properties%Elements(i)%Sigma_s(N%Group,N%Group))
        allocate(Properties%Elements(i)%Sigma_T(N%Group))
        allocate(Properties%Elements(i)%Sigma_r(N%Group))
        allocate(Properties%Elements(i)%Coordinates(this%N_Cell_Pointers(i),3))
        allocate(Properties%Elements(i)%Flux(N%Group,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Source(N%Group,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Source_Vector(this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Cell_Pointers(this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%D_Matrix(this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%A_Matrix(this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%B_Matrix(this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%K_Matrix(N%Group,this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
    end do

    do i = 1, N%Material

        read(1,*) ! Skip the line to align with the input file

        read(1,*) Properties%Materials(i)%Sigma_a 

        read(1,*) Properties%Materials(i)%Sigma_f

        do j = 1, N%Group

            read(1,*) sigma_s_vector

            Properties%Materials(i)%Sigma_s(j,:) = sigma_s_vector

            Properties%Materials(i)%Sigma_T(j) = Properties%Materials(i)%Sigma_a(j)+SUM(sigma_s_vector)

            Properties%Materials(i)%Sigma_r(j) = Properties%Materials(i)%Sigma_a(j)+SUM(sigma_s_vector)-sigma_s_vector(j)

        end do

        if (Properties%Adjoint == 1) then

            Properties%Materials(i)%Sigma_s = TRANSPOSE(Properties%Materials(i)%Sigma_s)

        end if

    end do

    do i = 1, N%Element

        Properties%Elements(i)%Material = this%Material_ID(i)

        Properties%Elements(i)%Cell_Type = this%Cell_Type(i)

        Properties%Elements(i)%Number_of_Nodes = this%N_Cell_Pointers(i)

        Properties%Elements(i)%Cell_Pointers = this%Cell_Pointers(i,:)

        if (this%Material_ID(i) > N%Material .or. this%Material_ID(i) < 1) then

            write(*,*) 'Invalid Material Index'

            stop

        end if

        Properties%Elements(i)%Sigma_a = Properties%Materials(this%Material_ID(i))%Sigma_a

        Properties%Elements(i)%Sigma_f = Properties%Materials(this%Material_ID(i))%Sigma_f

        Properties%Elements(i)%Sigma_s = Properties%Materials(this%Material_ID(i))%Sigma_s

        Properties%Elements(i)%Sigma_T = Properties%Materials(this%Material_ID(i))%Sigma_T

        Properties%Elements(i)%Sigma_r = Properties%Materials(this%Material_ID(i))%Sigma_r

        do j = 1, this%N_Cell_Pointers(i)
            Properties%Elements(i)%Coordinates(j,:) = this%Nodes(this%Cell_Pointers(i,j),:)
        end do

    end do

    N%N = this%N_nodes

    close(1)

end subroutine Read_Properties_Diffusion

end module m_Read_Properties_D