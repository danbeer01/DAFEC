module m_Construct_Matrix_3D
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
use m_Create_Hexahedral_Shape_Functions
use m_Create_Tetrahedral_Shape_Functions
use m_Create_Prismatic_Shape_Functions
use m_Create_Pyramidal_Shape_Functions
use m_VTK_Reader
use m_Gauss_Points

implicit none

contains

    subroutine Construct_Total_Matrix_3D(Properties, N, i, g_index, ang, mu, eta, xi)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        integer, intent(in)      :: i, g_index, ang
        real(kind=8), intent(in) :: mu, eta, xi

        real(kind = 8), dimension(:,:), allocatable   :: F_out

        allocate(F_out(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        call Construct_Streaming_Matrix(Properties, N, i, ang, mu, eta, xi)

        if (N%Degree > 1 .and. (Properties%Elements(1)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29)) then

            call Construct_F_out_Matrix_C(Properties, N, i, mu, eta, xi, F_out)

            call Construct_F_in_Matrix_C(Properties, N, i, ang, mu, eta, xi)

        else

            call Construct_F_out_Matrix(Properties, N, i, mu, eta, xi, F_out)

            call Construct_F_in_Matrix(Properties, N, i, ang, mu, eta, xi)

        end if

        Properties%Elements(i)%F_out_Matrix(ang,:,:) = F_out

        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix - Properties%Elements(i)%S_Matrix(ang,:,:) + F_out

    end subroutine Construct_Total_Matrix_3D

    subroutine Construct_Streaming_Matrix(Properties, N, i, ang, mu, eta, xi)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta, xi
        
        integer, intent(in)      :: i, ang

        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

            call Calculate_Hex_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta,xi)

        else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

            call Calculate_Tet_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta,xi)

        else if (Properties%Elements(i)%Cell_Type == 13) then

            call Calculate_Pris_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta,xi)

        else if (Properties%Elements(i)%Cell_Type == 14) then

            call Calculate_Pyr_Streaming_Matrix(Properties,N,i,Properties%Elements(i)%S_Matrix(ang,:,:),mu,eta,xi)

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Construct_Streaming_Matrix

    subroutine Construct_Mass_Matrix_3D(Properties, N, i)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N
        
        integer, intent(in) :: i

        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

            call Calculate_Hex_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

            call Calculate_Tet_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else if (Properties%Elements(i)%Cell_Type == 13) then

            call Calculate_Pris_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else if (Properties%Elements(i)%Cell_Type == 14) then

            call Calculate_Pyr_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Construct_Mass_Matrix_3D

    subroutine Construct_F_in_Matrix(Properties, N, i, ang, mu, eta, xi)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta, xi
        integer, intent(in)      :: i, ang
        
        real(kind = 8) :: Omega_n

        integer :: side_index

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = 0.0_8

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta + Properties%Elements(i)%Unit_Vectors(side_index,3)*xi)

            if (Omega_n < 0) then

                if (Properties%Elements(i)%Neighbours(side_index,1) == 0) then

                    if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Hex_Face(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tet_Face(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 13) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Pris_Face(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 14) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Pyr_Face(Properties,N,i,side_index)

                    end if

                else

                    if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Hex_Face_F_in(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tet_Face_F_in(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 13) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Pris_Face_F_in(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 14) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Pyr_Face_F_in(Properties,N,i,side_index)

                    end if

                end if

            end if

        end do

    end subroutine Construct_F_in_Matrix

    subroutine Construct_F_out_Matrix(Properties, N, i, mu, eta, xi, F_out)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta, xi
        integer, intent(in)      :: i

        real(kind = 8), dimension(:,:) :: F_out

        real(kind = 8) :: Omega_n

        integer :: side_index

        F_out = 0.0_8

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta + Properties%Elements(i)%Unit_Vectors(side_index,3)*xi)

            if (Omega_n > 0) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Hex_Face(Properties,N,i,side_index)

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Tet_Face(Properties,N,i,side_index)

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Pris_Face(Properties,N,i,side_index)

                else if (Properties%Elements(i)%Cell_Type == 14) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Pyr_Face(Properties,N,i,side_index)

                end if

            end if

        end do

    end subroutine Construct_F_out_Matrix

    subroutine Construct_F_in_Matrix_C(Properties, N, i, ang, mu, eta, xi)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta, xi
        integer, intent(in)      :: i, ang
        
        real(kind = 8) :: Omega_n

        integer :: side_index, gauss_index

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = 0.0_8

            do gauss_index = 1, (2*N%Degree+2)**2

                Omega_n = (Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,1)*mu + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,2)*eta + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,3)*xi)

                if (Omega_n < 0) then

                    if (Properties%Elements(i)%Neighbours(side_index,1) == 0) then

                        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                            call Integrate_Hex_Face_C(Properties,N,i,side_index,gauss_index,Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Omega_n)

                        end if

                    else

                        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                            call Integrate_Hex_Face_F_in_C(Properties,N,i,side_index,gauss_index,Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Omega_n)

                        end if

                    end if

                end if

            end do

        end do

    end subroutine Construct_F_in_Matrix_C

    subroutine Construct_F_out_Matrix_C(Properties, N, i, mu, eta, xi, F_out)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta, xi
        integer, intent(in)      :: i

        real(kind = 8), dimension(:,:) :: F_out

        real(kind = 8) :: Omega_n

        integer :: side_index, gauss_index

        F_out = 0.0_8

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            do gauss_index = 1, (2*N%Degree+2)**2

                Omega_n = (Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,1)*mu + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,2)*eta + Properties%Elements(i)%Gauss_Unit_Vectors(side_index,gauss_index,3)*xi)

                if (Omega_n > 0) then

                    if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                        call Integrate_Hex_Face_C(Properties,N,i,side_index,gauss_index,F_out,Omega_n)

                    end if

                end if

            end do

        end do

    end subroutine Construct_F_out_Matrix_C

    subroutine Calculate_Jacobian_3D(Properties, Degree, i)

        type(PropertiesType), intent(inout)  :: Properties

        integer, intent(in)                       :: Degree
        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w
        real(kind = 8)                            :: xi_val, eta_val, zeta_val
        integer, intent(in)                       :: i
        integer :: j, Num_Gauss_Points

        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

            Num_Gauss_Points = (Degree+1)**3

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_3D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)
                zeta_val = zeta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Hexahedral(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val, zeta_val)

                Properties%Elements(i)%Det_Jacobian(j) = Properties%Elements(i)%Jacobian(j,1,1)*(Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,2)) - Properties%Elements(i)%Jacobian(j,1,2)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,1)) + Properties%Elements(i)%Jacobian(j,1,3)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,2) - Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

            if (Degree == 1) then
                Num_Gauss_Points = 5
            else if (Degree == 2) then
                Num_Gauss_Points = 11
            end if

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_Tetrahedral_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)
                zeta_val = zeta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Tetrahedral(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val, zeta_val)

                Properties%Elements(i)%Det_Jacobian(j) = Properties%Elements(i)%Jacobian(j,1,1)*(Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,2)) - Properties%Elements(i)%Jacobian(j,1,2)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,1)) + Properties%Elements(i)%Jacobian(j,1,3)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,2) - Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else if (Properties%Elements(i)%Cell_Type == 13) then

            if (Degree == 1) then
                Num_Gauss_Points = 8
            else if (Degree == 2) then
                Num_Gauss_Points = 21
            end if

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_Prismatic_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)
                zeta_val = zeta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Prismatic(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val, zeta_val)

                Properties%Elements(i)%Det_Jacobian(j) = Properties%Elements(i)%Jacobian(j,1,1)*(Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,2)) - Properties%Elements(i)%Jacobian(j,1,2)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,1)) + Properties%Elements(i)%Jacobian(j,1,3)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,2) - Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else if (Properties%Elements(i)%Cell_Type == 14) then

            Num_Gauss_Points = 8

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,3,3), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_Pyramidal_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)
                zeta_val = zeta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Pyramidal(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val, zeta_val)

                Properties%Elements(i)%Det_Jacobian(j) = Properties%Elements(i)%Jacobian(j,1,1)*(Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,2)) - Properties%Elements(i)%Jacobian(j,1,2)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,3) - Properties%Elements(i)%Jacobian(j,2,3)*Properties%Elements(i)%Jacobian(j,3,1)) + Properties%Elements(i)%Jacobian(j,1,3)*(Properties%Elements(i)%Jacobian(j,2,1)*Properties%Elements(i)%Jacobian(j,3,2) - Properties%Elements(i)%Jacobian(j,2,2)*Properties%Elements(i)%Jacobian(j,3,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else

            write(*,*) 'Error: Cell type not supported'
            stop

        end if

    end subroutine Calculate_Jacobian_3D

    Function Jacobian_Hexahedral(this, Degree, Coordinates, eta, xi, zeta) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dz_deta, dx_dxi, dy_dxi, dz_dxi, dx_dzeta, dy_dzeta, dz_dzeta
    
        real(kind=8), intent(in) :: eta, xi, zeta
    
        real(kind=8), dimension(3,3) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dz_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8
        dz_deta = 0.0_8
        dx_dzeta = 0.0_8
        dy_dzeta = 0.0_8
        dz_dzeta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Hex_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Hex_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_deta = dz_deta + Generate_Hex_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dxi = dx_dxi + Generate_Hex_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Hex_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dxi = dz_dxi + Generate_Hex_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dzeta = dx_dzeta + Generate_Hex_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dzeta = dy_dzeta + Generate_Hex_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dzeta = dz_dzeta + Generate_Hex_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

        end do
    
        J(1,1) = dx_dxi
        J(2,1) = dy_dxi
        J(3,1) = dz_dxi
        J(1,2) = dx_deta
        J(2,2) = dy_deta
        J(3,2) = dz_deta
        J(1,3) = dx_dzeta
        J(2,3) = dy_dzeta
        J(3,3) = dz_dzeta
    
    end function Jacobian_Hexahedral

    Function Jacobian_Tetrahedral(this, Degree, Coordinates, eta, xi, zeta) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dz_deta, dx_dxi, dy_dxi, dz_dxi, dx_dzeta, dy_dzeta, dz_dzeta
    
        real(kind=8), intent(in) :: eta, xi, zeta
    
        real(kind=8), dimension(3,3) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dz_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8
        dz_deta = 0.0_8
        dx_dzeta = 0.0_8
        dy_dzeta = 0.0_8
        dz_dzeta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Tetrahedral_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Tetrahedral_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_deta = dz_deta + Generate_Tetrahedral_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dxi = dx_dxi + Generate_Tetrahedral_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Tetrahedral_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dxi = dz_dxi + Generate_Tetrahedral_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dzeta = dx_dzeta + Generate_Tetrahedral_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dzeta = dy_dzeta + Generate_Tetrahedral_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dzeta = dz_dzeta + Generate_Tetrahedral_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

        end do
    
        J(1,1) = dx_dxi
        J(2,1) = dy_dxi
        J(3,1) = dz_dxi
        J(1,2) = dx_deta
        J(2,2) = dy_deta
        J(3,2) = dz_deta
        J(1,3) = dx_dzeta
        J(2,3) = dy_dzeta
        J(3,3) = dz_dzeta
    
    end function Jacobian_Tetrahedral

    Function Jacobian_Prismatic(this, Degree, Coordinates, eta, xi, zeta) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dz_deta, dx_dxi, dy_dxi, dz_dxi, dx_dzeta, dy_dzeta, dz_dzeta
    
        real(kind=8), intent(in) :: eta, xi, zeta
    
        real(kind=8), dimension(3,3) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dz_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8
        dz_deta = 0.0_8
        dx_dzeta = 0.0_8
        dy_dzeta = 0.0_8
        dz_dzeta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Prismatic_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Prismatic_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_deta = dz_deta + Generate_Prismatic_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dxi = dx_dxi + Generate_Prismatic_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Prismatic_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dxi = dz_dxi + Generate_Prismatic_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dzeta = dx_dzeta + Generate_Prismatic_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dzeta = dy_dzeta + Generate_Prismatic_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dzeta = dz_dzeta + Generate_Prismatic_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

        end do

        J(1,1) = dx_dxi
        J(2,1) = dy_dxi
        J(3,1) = dz_dxi
        J(1,2) = dx_deta
        J(2,2) = dy_deta
        J(3,2) = dz_deta
        J(1,3) = dx_dzeta
        J(2,3) = dy_dzeta
        J(3,3) = dz_dzeta

    end function Jacobian_Prismatic

    Function Jacobian_Pyramidal(this, Degree, Coordinates, eta, xi, zeta) Result(J)

        type(PropertiesType), intent(in) :: this

        real(kind=8), dimension(:,:), intent(in) :: Coordinates
    
        real(kind=8) :: dx_deta, dy_deta, dz_deta, dx_dxi, dy_dxi, dz_dxi, dx_dzeta, dy_dzeta, dz_dzeta
    
        real(kind=8), intent(in) :: eta, xi, zeta
    
        real(kind=8), dimension(3,3) :: J

        integer, intent(in) :: Degree

        integer :: i

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dz_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8
        dz_deta = 0.0_8
        dx_dzeta = 0.0_8
        dy_dzeta = 0.0_8
        dz_dzeta = 0.0_8

        do i = 1, size(Coordinates,1)

            dx_deta = dx_deta + Generate_Pyramidal_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Pyramidal_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_deta = dz_deta + Generate_Pyramidal_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dxi = dx_dxi + Generate_Pyramidal_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Pyramidal_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dxi = dz_dxi + Generate_Pyramidal_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dzeta = dx_dzeta + Generate_Pyramidal_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dzeta = dy_dzeta + Generate_Pyramidal_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dzeta = dz_dzeta + Generate_Pyramidal_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

        end do

        J(1,1) = dx_dxi
        J(2,1) = dy_dxi
        J(3,1) = dz_dxi
        J(1,2) = dx_deta
        J(2,2) = dy_deta
        J(3,2) = dz_deta
        J(1,3) = dx_dzeta
        J(2,3) = dy_dzeta
        J(3,3) = dz_dzeta

    end function Jacobian_Pyramidal

    Function Inverse_Jacobian(J) Result(J_Inverse)

        real(kind=8), dimension(3,3) :: J
        
        real(kind=8), dimension(3,3) :: J_Inverse
    
        real(kind=8) :: Det_J
    
        Det_J = J(1,1)*(J(2,2)*J(3,3) - J(2,3)*J(3,2)) - J(1,2)*(J(2,1)*J(3,3) - J(2,3)*J(3,1)) + J(1,3)*(J(2,1)*J(3,2) - J(2,2)*J(3,1))
    
        J_Inverse(1,1) = (J(2,2)*J(3,3) - J(2,3)*J(3,2))/Det_J
        J_Inverse(1,2) = (J(1,3)*J(3,2) - J(1,2)*J(3,3))/Det_J
        J_Inverse(1,3) = (J(1,2)*J(2,3) - J(1,3)*J(2,2))/Det_J
        J_Inverse(2,1) = (J(2,3)*J(3,1) - J(2,1)*J(3,3))/Det_J
        J_Inverse(2,2) = (J(1,1)*J(3,3) - J(1,3)*J(3,1))/Det_J
        J_Inverse(2,3) = (J(1,3)*J(2,1) - J(1,1)*J(2,3))/Det_J
        J_Inverse(3,1) = (J(2,1)*J(3,2) - J(2,2)*J(3,1))/Det_J
        J_Inverse(3,2) = (J(1,2)*J(3,1) - J(1,1)*J(3,2))/Det_J
        J_Inverse(3,3) = (J(1,1)*J(2,2) - J(1,2)*J(2,1))/Det_J
    
      end function Inverse_Jacobian

end module m_Construct_Matrix_3D