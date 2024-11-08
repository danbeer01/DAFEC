module m_Construct_Matrix_1D
!
! Purpose:
! To construct the matrix of the linear system of equations for the neutron diffusion equation using the finite element method
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties
use m_Create_Shape_Functions
use m_VTK_Reader
use m_Gauss_Points

implicit none

contains

    subroutine Construct_Total_Matrix_1D(Properties, i, k, m, mu)

        type(PropertiesType), intent(inout)  :: Properties

        integer, intent(in)      :: i, k, m

        real(kind=8), intent(in) :: mu

        real(kind = 8), dimension(:,:), allocatable :: F_out

        integer :: side_index

        allocate(F_out(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        F_out = 0.0_8

        if (mu > 0.0_8) then

            F_out(2,2) = 1.0_8

        else if (mu < 0.0_8) then 

            F_out(1,1) = 1.0_8

        end if

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(m,:,:) = 0.0_8

            if (mu > 0.0_8) then

                if(side_index == 1) Properties%Elements(i)%Sides(side_index)%F_in_Matrix(m,1,2) = mu*1.0_8

            else if (mu < 0.0_8) then

                if(side_index == 2) Properties%Elements(i)%Sides(side_index)%F_in_Matrix(m,2,1) = abs(mu)*1.0_8

            end if

        end do

        if(Properties%g == 0) Properties%Elements(i)%K_Matrix(k,m,:,:) = -mu*Properties%Elements(i)%S_Matrix(1,:,:) + Properties%Elements(i)%Sigma_t(k)*Properties%Elements(i)%A_Matrix + ABS(mu)*F_out

    end subroutine Construct_Total_Matrix_1D

    subroutine Construct_Streaming_and_Mass_Matrix(Properties, N, i)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N
        
        integer, intent(in)      :: i

        call Calculate_Mass_Matrix(Properties, N, i, Properties%Elements(i)%A_Matrix)

        call Calculate_Streaming_Matrix(Properties, N, i, Properties%Elements(i)%S_Matrix(1,:,:))

    end subroutine Construct_Streaming_and_Mass_Matrix

    subroutine Calculate_Jacobian_1D(Properties, Degree, i)

        type(PropertiesType), intent(inout)  :: Properties

        integer, intent(in)                       :: Degree
        real(kind = 8), dimension(:), allocatable :: xi, w
        real(kind = 8)                            :: xi_val
        integer, intent(in)                       :: i
        integer :: j, Num_Gauss_Points

        if (Properties%g == 0) then
            Num_Gauss_Points = Degree+1
        else if (Properties%g == 1) then
            Num_Gauss_Points = Degree+2
        else if (Properties%g == 2) then
            Num_Gauss_Points = Degree+2
        end if

        allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,1,1), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,1,1), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        do j = 1, Num_Gauss_Points

            xi_val = xi(j)

            Properties%Elements(i)%Jacobian(j,:,:) = Jacobian(Properties, Degree, Properties%Elements(i)%Coordinates, xi_val)

            Properties%Elements(i)%Det_Jacobian(j) = Properties%Elements(i)%Jacobian(j,1,1)

            Properties%Elements(i)%Inverse_Jacobian(j,:,:) = 1.0_8/(Properties%Elements(i)%Jacobian(j,:,:))

        end do

        Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

    end subroutine Calculate_Jacobian_1D

    Function Jacobian(this, Degree, Coordinates, xi) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_dxi
    
        real(kind=8), intent(in) :: xi
    
        real(kind=8), dimension(1,1) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_dxi = dx_dxi + Generate_Shape_Functions_Derivative(this%Isoparametric_Coordinates(:,1), Degree, i, xi)*Coordinates(i,1)

        end do
    
        J(1,1) = dx_dxi
    
        end function Jacobian

end module m_Construct_Matrix_1D