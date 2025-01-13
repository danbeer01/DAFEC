module m_Construct_Matrix_2D
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
use m_Create_Quadrilateral_Shape_Functions
use m_Create_Triangular_Shape_Functions
use m_VTK_Reader
use m_Gauss_Points

implicit none

contains

    subroutine Construct_Total_Matrix_2D(Properties, N, i, g_index, ang, mu, eta)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        integer, intent(in)      :: i, g_index, ang
        real(kind=8), intent(in) :: mu, eta

        real(kind = 8), dimension(:,:), allocatable   :: F_out

        allocate(F_out(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        call Construct_Streaming_Matrix(Properties, N, i, ang, mu, eta)

        if (N%Degree == 1) then

            call Construct_F_out_Matrix(Properties, N, i, mu, eta, F_out)

            call Construct_F_in_Matrix(Properties, N, i, ang, mu, eta)

        else

            call Construct_F_out_Matrix_C(Properties, N, i, mu, eta, F_out)

            call Construct_F_in_Matrix_C(Properties, N, i, ang, mu, eta)

        end if

        Properties%Elements(i)%F_out_Matrix(ang,:,:) = F_out

        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = -Properties%Elements(i)%S_Matrix(ang,:,:) + Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix + F_out

    end subroutine Construct_Total_Matrix_2D

    subroutine Construct_Streaming_Matrix(Properties, N, i, ang, mu, eta)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        
        integer, intent(in)      :: i, ang

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

            call Calculate_Tri_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta)

        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

            call Calculate_Quad_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta)

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Construct_Streaming_Matrix

    subroutine Construct_Mass_Matrix_2D(Properties, N, i)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N
        
        integer, intent(in) :: i

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

            call Calculate_Tri_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

            call Calculate_Quad_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Construct_Mass_Matrix_2D

    subroutine Construct_F_in_Matrix(Properties, N, i, ang, mu, eta)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i, ang
        
        real(kind = 8) :: Omega_n

        integer :: side_index

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = 0.0_8

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta)

            if (Omega_n < 0) then

                if (Properties%Elements(i)%Neighbours(side_index,1) == 0) then

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tri_Side(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Quad_Side(Properties,N,i,side_index)

                    end if

                else

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tri_Side_F_in(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Quad_Side_F_in(Properties,N,i,side_index)

                    end if

                end if

            end if

        end do

    end subroutine Construct_F_in_Matrix

    subroutine Construct_F_out_Matrix(Properties, N, i, mu, eta, F_out)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i

        real(kind = 8), dimension(:,:) :: F_out

        real(kind = 8) :: Omega_n

        integer :: side_index

        F_out = 0.0_8

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta)

            if (Omega_n > 0) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Tri_Side(Properties,N,i,side_index)

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Quad_Side(Properties,N,i,side_index)

                end if

            end if

        end do

    end subroutine Construct_F_out_Matrix

    subroutine Construct_F_in_Matrix_C(Properties, N, i, ang, mu, eta)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i, ang
        
        real(kind = 8) :: Omega_n

        integer :: side_index, gauss_index

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = 0.0_8

            do gauss_index = 1, 2*N%Degree + 2

                Omega_n = (Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,1)*mu + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,2)*eta)

                if (Omega_n < 0) then

                    if (Properties%Elements(i)%Neighbours(side_index,1) == 0) then

                        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                            print *, 'Triangular cell type not supported'
                            stop

                        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 23 .or. Properties%Elements(i)%Cell_Type == 28) then

                            call Integrate_Quad_Side_C(Properties,N,i,side_index,gauss_index,Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Omega_n)

                        end if

                    else

                        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                            print *, 'Triangular cell type not supported'
                            stop

                        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 23 .or. Properties%Elements(i)%Cell_Type == 28) then

                            call Integrate_Quad_Side_F_in_C(Properties,N,i,side_index,gauss_index,Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Omega_n)

                        end if

                    end if

                end if

            end do

        end do

    end subroutine Construct_F_in_Matrix_C

    subroutine Construct_F_out_Matrix_C(Properties, N, i, mu, eta, F_out)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i

        real(kind = 8), dimension(:,:) :: F_out

        real(kind = 8) :: Omega_n, r

        integer :: side_index, gauss_index

        F_out = 0.0_8

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1) + Properties%Elements(i)%Coordinates(3,1) + Properties%Elements(i)%Coordinates(4,1))/4.0_8

            do gauss_index = 1, 2*N%Degree + 2

                Omega_n = (Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,1)*mu + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,2)*eta)

                if (Omega_n > 0) then

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                        print *, 'Triangular cell type not supported'
                        stop

                    else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 23 .or. Properties%Elements(i)%Cell_Type == 28) then

                        call Integrate_Quad_Side_C(Properties,N,i,side_index,gauss_index,F_out,Omega_n)

                    end if

                end if

            end do

        end do

        if (Properties%g == 1) F_out = F_out*r

    end subroutine Construct_F_out_Matrix_C

    subroutine Calculate_Jacobian_2D(Properties, Degree, i)

        type(PropertiesType), intent(inout)  :: Properties

        integer, intent(in)                       :: Degree
        real(kind = 8), dimension(:), allocatable :: xi, w, eta
        real(kind = 8)                            :: xi_val, eta_val
        integer, intent(in)                       :: i
        integer :: j, Num_Gauss_Points

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type  == 22) then

            if(size(Properties%Elements(i)%Coordinates,1) == 3) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Triangular(Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val)

                Properties%Elements(i)%Det_Jacobian(j) = (Properties%Elements(i)%Jacobian(j,1,1)*Properties%Elements(i)%Jacobian(j,2,2)) - (Properties%Elements(i)%Jacobian(j,1,2)*Properties%Elements(i)%Jacobian(j,2,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = 0.5_8*SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Quadrilateral(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val)

                Properties%Elements(i)%Det_Jacobian(j) = (Properties%Elements(i)%Jacobian(j,1,1)*Properties%Elements(i)%Jacobian(j,2,2)) - (Properties%Elements(i)%Jacobian(j,1,2)*Properties%Elements(i)%Jacobian(j,2,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Calculate_Jacobian_2D

    Function Jacobian_Quadrilateral(this, Degree, Coordinates, eta, xi) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dx_dxi, dy_dxi
    
        real(kind=8), intent(in) :: eta, xi
    
        real(kind=8), dimension(2,2) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi)*Coordinates(i,2)

            dx_dxi = dx_dxi + Generate_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi)*Coordinates(i,2)

        end do
    
        J(1,1) = dx_dxi
        J(1,2) = dy_dxi
        J(2,1) = dx_deta
        J(2,2) = dy_deta
    
    end function Jacobian_Quadrilateral

    Function Jacobian_Triangular(Degree, Coordinates, eta, xi) Result(J)
    
        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dx_dxi, dy_dxi
    
        real(kind=8), intent(in) :: eta, xi
    
        real(kind=8), dimension(2,2) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Tri_Shape_Functions_Derivative_eta(Degree, i, eta, xi)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Tri_Shape_Functions_Derivative_eta(Degree, i, eta, xi)*Coordinates(i,2)

            dx_dxi = dx_dxi + Generate_Tri_Shape_Functions_Derivative_xi(Degree, i, eta, xi)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Tri_Shape_Functions_Derivative_xi(Degree, i, eta, xi)*Coordinates(i,2)

        end do
    
        J(1,1) = dx_dxi
        J(1,2) = dy_dxi
        J(2,1) = dx_deta
        J(2,2) = dy_deta
    
    end function Jacobian_Triangular

    Function Inverse_Jacobian(J) Result(J_Inverse)

        real(kind=8), dimension(2,2) :: J
        
        real(kind=8), dimension(2,2) :: J_Inverse
    
        real(kind=8) :: Det_J
    
        Det_J = J(1,1)*J(2,2) - J(1,2)*J(2,1)
    
        J_Inverse(1,1) = J(2,2)/Det_J
        J_Inverse(1,2) = -J(1,2)/Det_J
        J_Inverse(2,1) = -J(2,1)/Det_J
        J_Inverse(2,2) = J(1,1)/Det_J
    
    end function Inverse_Jacobian

end module m_Construct_Matrix_2D