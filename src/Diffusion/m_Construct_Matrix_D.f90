module m_Construct_Matrix_D
!
! Purpose:
! To construct the matrix of the linear system of equations for the neutron diffusion equation using the finite element method
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Create_Hexahedral_Shape_Functions_D
use m_Create_Tetrahedral_Shape_Functions_D
use m_Create_Prismatic_Shape_Functions_D
use m_Create_Pyramidal_Shape_Functions_D
use m_Create_Quadrilateral_Shape_Functions_D
use m_Create_Triangular_Shape_Functions_D
use m_Create_Shape_Functions_D
use m_VTK_Reader
use m_Boundary_Conditions_D
use m_Gauss_Points

implicit none

contains

    subroutine Construct_Matrix(Properties, N, Periodic_Pairs)

        type(PropertiesTypeD), intent(inout)  :: Properties
        type(NTypeD), intent(in)              :: N

        integer :: i, k

        integer, dimension(:,:), allocatable :: Periodic_Pairs

        if(N%D == 2) call Calculate_Isoparametric_Quadrilateral_Coordinates(N%Degree, Properties)
        if(N%D == 3) call Calculate_Isoparametric_Hexahedral_Coordinates(N%Degree, Properties)

        do i = 1, N%Element

            if(N%D == 1) call Calculate_Jacobian_1D(Properties, i)
            if(N%D == 2) call Calculate_Jacobian_2D(Properties, N%Degree, i)
            if(N%D == 3) call Calculate_Jacobian_3D(Properties, N%Degree, i)

            if (Properties%Elements(i)%Cell_Type == 3 .or. Properties%Elements(i)%Cell_Type == 21) then

                call Calculate_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                call Calculate_Tri_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Tri_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Tri_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                call Calculate_Quad_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Quad_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Quad_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                call Calculate_Hex_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Hex_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Hex_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

                call Calculate_Tet_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Tet_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Tet_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 13) then

                call Calculate_Pris_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Pris_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Pris_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else if (Properties%Elements(i)%Cell_Type == 14) then

                call Calculate_Pyr_D_Matrix(Properties,N,i,Properties%Elements(i)%D_Matrix)

                call Calculate_Pyr_Mass_Matrix(Properties,N,i,Properties%Elements(i)%A_Matrix)

                call Calculate_Pyr_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

            else

                write(*,*) 'Error: Cell type not supported'
                stop

            end if

            if(N%D == 1) call Calculate_B_Matrix_1D(Properties,i,Properties%Elements(i)%B_Matrix)
            if(N%D == 2) call Calculate_B_Matrix_2D(Properties,N,i,Properties%Elements(i)%B_Matrix)
            if(N%D == 3) call Calculate_B_Matrix_3D(Properties,N,i,Properties%Elements(i)%B_Matrix)

        end do

        do i = 1, N%Element

            do k = 1, N%Group

                Properties%Elements(i)%K_Matrix(k,:,:) = Properties%Elements(i)%Sigma_r(k)*Properties%Elements(i)%A_Matrix + (1.0_8/(3.0_8*Properties%Elements(i)%Sigma_t(k)))*Properties%Elements(i)%D_Matrix + Properties%Elements(i)%B_Matrix

            end do

        end do

        if((N%D == 2) .and. (Properties%LBC == 4 .or. Properties%RBC == 4 .or. Properties%TBC == 4 .or. Properties%BBC == 4)) call Periodic_Boundary_2D(N,Properties,Periodic_Pairs)
        if((N%D == 3) .and. (Properties%LBC == 4 .or. Properties%RBC == 4 .or. Properties%TBC == 4 .or. Properties%BBC == 4 .or. Properties%FBC == 4.or. Properties%BBC == 4)) call Periodic_Boundary_3D(N,Properties,Periodic_Pairs)

    end subroutine Construct_Matrix

    subroutine Calculate_Jacobian_1D(Properties, i)

        type(PropertiesTypeD), intent(inout) :: Properties

        integer, intent(in) :: i
        integer             :: Num_Gauss_Points

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points - 1

        allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,1,1), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,1,1), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

        Properties%Elements(i)%Volume = ABS(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))
        
    end subroutine Calculate_Jacobian_1D

    subroutine Calculate_Jacobian_2D(Properties, Degree, i)

        type(PropertiesTypeD), intent(inout)  :: Properties

        integer, intent(in)                       :: Degree
        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        real(kind = 8)                            :: xi_val, eta_val
        integer, intent(in)                       :: i
        integer :: j, Num_Gauss_Points

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Triangular(Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val)

                Properties%Elements(i)%Det_Jacobian(j) = (Properties%Elements(i)%Jacobian(j,1,1)*Properties%Elements(i)%Jacobian(j,2,2)) - (Properties%Elements(i)%Jacobian(j,1,2)*Properties%Elements(i)%Jacobian(j,2,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian_2D(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = 0.5_8*SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

            allocate(Properties%Elements(i)%Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Inverse_Jacobian(Num_Gauss_Points,2,2), Properties%Elements(i)%Det_Jacobian(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

            call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            do j = 1, Num_Gauss_Points

                xi_val = xi(j)
                eta_val = eta(j)

                Properties%Elements(i)%Jacobian(j,:,:) = Jacobian_Quadrilateral(Properties, Degree, Properties%Elements(i)%Coordinates, eta_val, xi_val)

                Properties%Elements(i)%Det_Jacobian(j) = (Properties%Elements(i)%Jacobian(j,1,1)*Properties%Elements(i)%Jacobian(j,2,2)) - (Properties%Elements(i)%Jacobian(j,1,2)*Properties%Elements(i)%Jacobian(j,2,1))

                Properties%Elements(i)%Inverse_Jacobian(j,:,:) = Inverse_Jacobian_2D(Properties%Elements(i)%Jacobian(j,:,:))

            end do

            Properties%Elements(i)%Volume = SUM(w*ABS(Properties%Elements(i)%Det_Jacobian))

        end if
        
    end subroutine Calculate_Jacobian_2D

    Function Jacobian_Quadrilateral(this, Degree, Coordinates, eta, xi) Result(J)

        type(PropertiesTypeD), intent(in) :: this

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

            dx_deta = dx_deta + Generate_Quad_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Quad_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi)*Coordinates(i,2)

            dx_dxi = dx_dxi + Generate_Quad_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Quad_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi)*Coordinates(i,2)

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

    subroutine Calculate_Jacobian_3D(Properties, Degree, i)

        type(PropertiesTypeD), intent(inout)  :: Properties

        integer, intent(in)                       :: Degree
        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w
        real(kind = 8)                            :: xi_val, eta_val, zeta_val
        integer, intent(in)                       :: i
        integer :: j, Num_Gauss_Points

        if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

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

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

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

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

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

            Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

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

        type(PropertiesTypeD), intent(in) :: this

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

            dx_deta = dx_deta + Generate_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_deta = dy_deta + Generate_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_deta = dz_deta + Generate_Shape_Functions_Derivative_eta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dxi = dx_dxi + Generate_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dxi = dy_dxi + Generate_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dxi = dz_dxi + Generate_Shape_Functions_Derivative_xi(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

            dx_dzeta = dx_dzeta + Generate_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,1)

            dy_dzeta = dy_dzeta + Generate_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,2)

            dz_dzeta = dz_dzeta + Generate_Shape_Functions_Derivative_zeta(this, Degree, i, eta, xi, zeta)*Coordinates(i,3)

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

        type(PropertiesTypeD), intent(in) :: this

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

        type(PropertiesTypeD), intent(in) :: this

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

        type(PropertiesTypeD), intent(in) :: this

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

      Function Inverse_Jacobian_2D(J) Result(J_Inverse)

        real(kind=8), dimension(2,2) :: J
        
        real(kind=8), dimension(2,2) :: J_Inverse
    
        real(kind=8) :: Det_J
    
        Det_J = J(1,1)*J(2,2) - J(1,2)*J(2,1)
    
        J_Inverse(1,1) = J(2,2)/Det_J
        J_Inverse(1,2) = -J(1,2)/Det_J
        J_Inverse(2,1) = -J(2,1)/Det_J
        J_Inverse(2,2) = J(1,1)/Det_J
    
      end function Inverse_Jacobian_2D

end module