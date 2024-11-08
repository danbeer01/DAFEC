module m_Read_Properties
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

type :: NType
    integer :: N
    integer :: D
    integer :: Group
    integer :: Degree
    integer :: Element
    integer :: Material
    integer :: Angle
    integer :: Ordinates
    integer :: Anisotropy
end type NType

type :: MaterialType
    real(kind=8), dimension(:), allocatable     :: Sigma_a
    real(kind=8), dimension(:), allocatable     :: Sigma_f
    real(kind=8), dimension(:,:,:), allocatable :: Sigma_s
    real(kind=8), dimension(:), allocatable     :: Sigma_t
end type MaterialType

type :: SideType
    real(kind=8), dimension(:,:,:), allocatable :: F_in_Matrix
    real(kind=8), dimension(:,:,:), allocatable :: Boundary
    integer :: Nodes_on_Side
    integer, dimension(:), allocatable :: Neighbour_Nodes
end type SideType

type :: ElementType
    integer :: Material
    integer, dimension(:), allocatable :: Cell_Pointers
    integer :: Number_of_Nodes
    integer :: Number_of_Sides
    integer :: Cell_Type
    real(kind=8) :: Volume
    real(kind=8), dimension(:,:), allocatable     :: Coordinates
    real(kind=8), dimension(:), allocatable       :: Sigma_a
    real(kind=8), dimension(:), allocatable       :: Sigma_f
    real(kind=8), dimension(:,:,:), allocatable   :: Sigma_s
    real(kind=8), dimension(:), allocatable       :: Sigma_t
    real(kind=8), dimension(:,:,:), allocatable   :: Jacobian
    real(kind=8), dimension(:,:,:), allocatable   :: Inverse_Jacobian
    real(kind=8), dimension(:), allocatable       :: Det_Jacobian
    real(kind=8), dimension(:,:,:,:), allocatable :: K_Matrix
    real(kind=8), dimension(:,:), allocatable     :: A_Matrix
    real(kind=8), dimension(:,:,:), allocatable   :: S_Matrix
    real(kind=8), dimension(:,:,:), allocatable   :: Source
    real(kind=8), dimension(:,:,:), allocatable   :: Sweep_Source
    real(kind=8), dimension(:,:,:), allocatable   :: Total_Source
    real(kind=8), dimension(:), allocatable       :: Source_Vector
    real(kind=8), dimension(:,:,:), allocatable   :: Flux
    real(kind=8), dimension(:,:), allocatable     :: Scalar_Flux
    real(kind=8), dimension(:,:,:,:), allocatable :: Legendre_Flux
    real(kind=8), dimension(:,:), allocatable     :: Unit_Vectors
    integer, dimension(:,:), allocatable          :: Neighbours
    integer, dimension(:,:), allocatable          :: Side_Nodes
    type(SideType), dimension(:), allocatable     :: Sides
end type ElementType

type :: PropertiesType
    integer :: LBC,RBC,TBC,BBC,FBC,BaBC
    integer :: Adjoint
    integer :: Case
    integer :: g
    real(kind = 8) :: alpha, Q_s
    real(kind = 8), dimension(:), allocatable     :: Length
    real(kind = 8), dimension(:), allocatable     :: Chi
    real(kind = 8), dimension(:,:), allocatable   :: Isoparametric_Coordinates
    type(ElementType), dimension(:), allocatable  :: Elements
    type(MaterialType), dimension(:), allocatable :: Materials
end type PropertiesType

contains

subroutine Read_Properties(Properties, N, this)
    implicit none

    class(Mesh)                          :: this
    class(PropertiesType), intent(inout) :: Properties
    class(NType), intent(inout) :: N
    
    integer :: i, j, l, Number_of_Sides
    real(kind=8), dimension(:), allocatable :: sigma_s_vector

    read(1,*)
    read(1,*) Properties%Case
    read(1,*) N%Material
    read(1,*) N%Group
    read(1,*) N%Angle
    read(1,*) N%Anisotropy
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
    N%D = this%Dimension
    if(N%D == 1) then
        if(Properties%g == 0 .or. Properties%g == 2) N%Ordinates = N%Angle
        if(Properties%g == 1) N%Ordinates = (N%Angle+2)*N%Angle/4
    else if(N%D == 2) then
        N%Ordinates = (N%Angle+2)*N%Angle/2
    else if(N%D == 3) then
        N%Ordinates = (N%Angle+2)*N%Angle
    end if
    
    allocate(Properties%Length(N%D))
    allocate(Properties%Chi(N%Group))
    
    read(1,*) ! Skip a line to align with input file
    read(1,*) Properties%Chi

    if ((SUM(Properties%Chi) - 1.0_8) > 1e-4) then
        print *, Properties%Chi
        write(*,*) 'Chi values must sum to 1'
        stop
    end if

    N%Element = this%N_cells

    do i = 1, N%D
        Properties%Length(i) = MAXVAL(this%Nodes(:,i)) - MINVAL(this%Nodes(:,i))
    end do
    
    allocate(Properties%Elements(N%Element))
    allocate(Properties%Materials(N%Material))
    allocate(sigma_s_vector(N%Group))

    do i = 1, N%Material
        allocate(Properties%Materials(i)%Sigma_a(N%Group))
        allocate(Properties%Materials(i)%Sigma_f(N%Group))
        allocate(Properties%Materials(i)%Sigma_s(0:N%Anisotropy,N%Group,N%Group))
        allocate(Properties%Materials(i)%Sigma_T(N%Group))
    end do

    do i = 1, N%Element
        allocate(Properties%Elements(i)%Sigma_a(N%Group))
        allocate(Properties%Elements(i)%Sigma_f(N%Group))
        allocate(Properties%Elements(i)%Sigma_s(0:N%Anisotropy,N%Group,N%Group))
        allocate(Properties%Elements(i)%Sigma_T(N%Group))
        allocate(Properties%Elements(i)%Coordinates(this%N_Cell_Pointers(i),3))
        allocate(Properties%Elements(i)%Source(N%Group,N%Ordinates,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Sweep_Source(N%Group,N%Ordinates,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Total_Source(N%Group,N%Ordinates,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Source_Vector(this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Scalar_Flux(N%Group,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Flux(N%Group,N%Ordinates,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Legendre_Flux(N%Group,0:N%Anisotropy,0:N%Anisotropy,this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%S_Matrix(N%Ordinates,this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%A_Matrix(this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        if (this%Cell_Type(i) == 3 .or. this%Cell_Type(i) == 21) then

            Number_of_Sides = 2
            Properties%Elements(i)%Number_of_Sides = 2
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))

        else if (this%Cell_Type(i) == 9 .or. this%Cell_Type(i) == 28) then

            Number_of_Sides = 4
            Properties%Elements(i)%Number_of_Sides = 4
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,N%Degree+1))

        else if (this%Cell_Type(i) == 5 .or. this%Cell_Type(i) == 22) then

            Number_of_Sides = 3
            Properties%Elements(i)%Number_of_Sides = 3
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,N%Degree+1))

        else if (this%Cell_Type(i) == 12 .or. this%Cell_Type(i) == 29) then

            Number_of_Sides = 6
            Properties%Elements(i)%Number_of_Sides = Number_of_Sides
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,(N%Degree+1)**2))
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            Properties%Elements(i)%Sides(:)%Nodes_on_Side = 4

        else if (this%Cell_Type(i) == 10 .or. this%Cell_Type(i) == 24) then

            Number_of_Sides = 4
            Properties%Elements(i)%Number_of_Sides = Number_of_Sides
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,(N%Degree+1)*(N%Degree+2)/2))
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            Properties%Elements(i)%Sides(:)%Nodes_on_Side = 3

        else if (this%Cell_Type(i) == 13) then

            Number_of_Sides = 5
            Properties%Elements(i)%Number_of_Sides = Number_of_Sides
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,(N%Degree)*4))
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            Properties%Elements(i)%Sides(1)%Nodes_on_Side = 4
            Properties%Elements(i)%Sides(2)%Nodes_on_Side = 4
            Properties%Elements(i)%Sides(3)%Nodes_on_Side = 4
            Properties%Elements(i)%Sides(4)%Nodes_on_Side = 3
            Properties%Elements(i)%Sides(5)%Nodes_on_Side = 3

        else if (this%Cell_Type(i) == 14) then

            Number_of_Sides = 5
            Properties%Elements(i)%Number_of_Sides = Number_of_Sides
            allocate(Properties%Elements(i)%Side_Nodes(Number_of_Sides,(N%Degree+1)**2))
            allocate(Properties%Elements(i)%Sides(Number_of_Sides))
            Properties%Elements(i)%Sides(1)%Nodes_on_Side = 4
            Properties%Elements(i)%Sides(2)%Nodes_on_Side = 3
            Properties%Elements(i)%Sides(3)%Nodes_on_Side = 3
            Properties%Elements(i)%Sides(4)%Nodes_on_Side = 3
            Properties%Elements(i)%Sides(5)%Nodes_on_Side = 3

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if
        allocate(Properties%Elements(i)%K_Matrix(N%Group,N%Ordinates,this%N_Cell_Pointers(i),this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Cell_Pointers(this%N_Cell_Pointers(i)))
        allocate(Properties%Elements(i)%Neighbours(Number_of_Sides,2))
    end do

    do i = 1, N%Material

        read(1,*) ! Skip the line to align with the input file

        read(1,*) Properties%Materials(i)%Sigma_a 

        read(1,*) Properties%Materials(i)%Sigma_f

        do l = 0, N%Anisotropy

            do j = 1, N%Group

                read(1,*) sigma_s_vector

                Properties%Materials(i)%Sigma_s(l,j,:) = sigma_s_vector

                if (l == 0) then

                    Properties%Materials(i)%Sigma_T(j) = Properties%Materials(i)%Sigma_a(j)+SUM(sigma_s_vector)

                end if

            end do

            if (Properties%Adjoint == 1) then

                Properties%Materials(i)%Sigma_s(l,:,:) = TRANSPOSE(Properties%Materials(i)%Sigma_s(l,:,:))

            end if

        end do

    end do

    do i = 1, N%Element

        Properties%Elements(i)%Material = this%Material_ID(i)

        Properties%Elements(i)%Cell_Type = this%Cell_Type(i)

        Properties%Elements(i)%Cell_Pointers = this%Cell_Pointers(i,:)

        Properties%Elements(i)%Number_of_Nodes = this%N_Cell_Pointers(i)

        if (this%Material_ID(i) > N%Material .or. this%Material_ID(i) < 1) then

            write(*,*) 'Invalid Material Index'

            stop

        end if

        Properties%Elements(i)%Sigma_a = Properties%Materials(this%Material_ID(i))%Sigma_a

        Properties%Elements(i)%Sigma_f = Properties%Materials(this%Material_ID(i))%Sigma_f

        Properties%Elements(i)%Sigma_s = Properties%Materials(this%Material_ID(i))%Sigma_s

        Properties%Elements(i)%Sigma_T = Properties%Materials(this%Material_ID(i))%Sigma_T

        do j = 1, this%N_Cell_Pointers(i)
            Properties%Elements(i)%Coordinates(j,:) = this%Nodes(this%Cell_Pointers(i,j),:)
        end do

    end do

    N%N = this%N_nodes

    close(1)

end subroutine Read_Properties

end module m_Read_Properties