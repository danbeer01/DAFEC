module m_Create_Triangular_Shape_Functions_D
!
! Purpose:
! To create the shape functions for a triangular element
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Gauss_Points

implicit none

contains

    Function Generate_Triangular_Shape_Functions(Degree, a, eta, xi) Result(Value)

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: eta, xi

        allocate(Shape_Functions(((Degree+1)*(Degree+2))/2))

        if (Degree == 1) then

            Shape_Functions(1) = 1.0_8 - xi - eta
            Shape_Functions(2) = xi
            Shape_Functions(3) = eta

        else if (Degree == 2) then

            Shape_Functions(1) = 2.0_8*(1.0_8 - xi - eta)*(0.5_8 - xi - eta)
            Shape_Functions(2) = 2.0_8*xi*(xi - 0.5_8)
            Shape_Functions(3) = 2.0_8*eta*(eta - 0.5_8)
            Shape_Functions(4) = 4.0_8*(1.0_8 - xi - eta)*xi
            Shape_Functions(5) = 4.0_8*xi*eta
            Shape_Functions(6) = 4.0_8*(1.0_8 - xi - eta)*eta

        else if (Degree == 3) then

            Shape_Functions(1) = 4.5_8*(1.0_8 - xi - eta)*(1.0_8/3.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions(2) = 4.5_8*xi*(xi - 1.0_8/3.0_8)*(xi - 2.0_8/3.0_8)
            Shape_Functions(3) = 4.5_8*eta*(eta - 1.0_8/3.0_8)*(eta - 2.0_8/3.0_8)
            Shape_Functions(4) = 13.5_8*(1.0_8 - xi - eta)*xi*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions(5) = 13.5_8*(1.0_8 - xi - eta)*xi*(xi - 1.0_8/3.0_8)
            Shape_Functions(6) = 13.5_8*xi*eta*(xi - 1.0_8/3.0_8)
            Shape_Functions(7) = -13.5_8*xi*eta*(1.0_8/3.0_8 - eta)
            Shape_Functions(8) = 13.5_8*(1.0_8 - xi - eta)*eta*(eta - 1.0_8/3.0_8)
            Shape_Functions(9) = 13.5_8*(1.0_8 - xi - eta)*eta*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions(10) = 27.0_8*(1.0_8 - xi - eta)*xi*eta

        else

            print*, 'Triangular elements only supported up to third order'
            stop

        end if

        Value = Shape_Functions(a)

    end Function Generate_Triangular_Shape_Functions

    Function Generate_Tri_Shape_Functions_Derivative_xi(Degree, a, eta, xi) Result(Value)

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta

        allocate(Shape_Functions_Derivative(((Degree+1)*(Degree+2))/2))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = -1.0_8
            Shape_Functions_Derivative(2) = 1.0_8
            Shape_Functions_Derivative(3) = 0.0_8

        else if (Degree == 2) then

            Shape_Functions_Derivative(1) = -3.0_8 + 4.0_8*xi + 4.0_8*eta
            Shape_Functions_Derivative(2) = 4.0_8*xi - 1.0_8
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 4.0_8 - 8.0_8*xi - 4.0_8*eta
            Shape_Functions_Derivative(5) = 4.0_8*eta
            Shape_Functions_Derivative(6) = -4.0_8*eta

        else if (Degree == 3) then

            Shape_Functions_Derivative(1) = -4.5_8*(1.0_8 - xi - eta)*(1.0_8/3.0_8 - xi - eta) - 4.5_8*(1.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta) - 4.5_8*(1.0_8/3.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions_Derivative(2) = 4.5_8*xi*(xi - 1.0_8/3.0_8) + 4.5_8*xi*(xi - 2.0_8/3.0_8) + 4.5_8*(xi - 1.0_8/3.0_8)*(xi - 2.0_8/3.0_8)
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 13.5_8*(1.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta) - 13.5_8*(1.0_8 - xi - eta)*xi - 13.5_8*xi*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions_Derivative(5) = 13.5_8*(1.0_8 - xi - eta)*(xi - 1.0_8/3.0_8) + 13.5_8*(1.0_8 - xi - eta)*xi - 13.5_8*xi*(xi - 1.0_8/3.0_8)
            Shape_Functions_Derivative(6) = 13.5_8*eta*(xi - 1.0_8/3.0_8) + 13.5_8*eta*xi
            Shape_Functions_Derivative(7) = -13.5_8*eta*(1.0_8/3.0_8 - eta)
            Shape_Functions_Derivative(8) = -13.5_8*eta*(eta - 1.0_8/3.0_8)
            Shape_Functions_Derivative(9) = -13.5_8*eta*(2.0_8/3.0_8 - xi - eta) - 13.5_8*eta*(1.0_8 - xi - eta)
            Shape_Functions_Derivative(10) = 27.0_8*(1.0_8 - xi - eta)*eta - 27.0_8*xi*eta

        else

            print*, 'Triangular elements only supported up to third order'
            stop

        end if

        Value = Shape_Functions_Derivative(a)

    end Function Generate_Tri_Shape_Functions_Derivative_xi

    Function Generate_Tri_Shape_Functions_Derivative_eta(Degree, a, eta, xi) Result(Value)

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta

        allocate(Shape_Functions_Derivative(((Degree+1)*(Degree+2))/2))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = -1.0_8
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 1.0_8

        else if (Degree == 2) then

            Shape_Functions_Derivative(1) = -3.0_8 + 4.0_8*xi + 4.0_8*eta
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 4.0_8*eta - 1.0_8
            Shape_Functions_Derivative(4) = -4.0_8*xi
            Shape_Functions_Derivative(5) = 4.0_8*xi
            Shape_Functions_Derivative(6) = 4.0_8 - 4.0_8*xi - 8.0_8*eta

        else if (Degree == 3) then

            Shape_Functions_Derivative(1) = -4.5_8*(1.0_8 - xi - eta)*(1.0_8/3.0_8 - xi - eta) - 4.5_8*(1.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta) - 4.5_8*(1.0_8/3.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 4.5_8*eta*(eta - 1.0_8/3.0_8) + 4.5_8*eta*(eta - 2.0_8/3.0_8) + 4.5_8*(eta - 1.0_8/3.0_8)*(eta - 2.0_8/3.0_8)
            Shape_Functions_Derivative(4) = -13.5_8*xi*(2.0_8/3.0_8 - xi - eta) - 13.5_8*xi*(1.0_8 - xi - eta)
            Shape_Functions_Derivative(5) = -13.5_8*xi*(xi - 1.0_8/3.0_8)
            Shape_Functions_Derivative(6) = 13.5_8*xi*(xi - 1.0_8/3.0_8)
            Shape_Functions_Derivative(7) = -13.5_8*xi*(1.0_8/3.0_8 - eta) + 13.5_8*xi*eta
            Shape_Functions_Derivative(8) = 13.5_8*(1.0_8 - xi - eta)*(eta - 1.0_8/3.0_8) + 13.5_8*(1.0_8 - xi - eta)*eta - 13.5_8*eta*(eta - 1.0_8/3.0_8)
            Shape_Functions_Derivative(9) = 13.5_8*(1.0_8 - xi - eta)*(2.0_8/3.0_8 - xi - eta) - 13.5_8*(1.0_8 - xi - eta)*eta - 13.5_8*eta*(2.0_8/3.0_8 - xi - eta)
            Shape_Functions_Derivative(10) = 27.0_8*(1.0_8 - xi - eta)*xi - 27.0_8*xi*eta

        else

            print*, 'Triangular elements only supported up to third order'
            stop

        end if

        Value = Shape_Functions_Derivative(a)

    end Function Generate_Tri_Shape_Functions_Derivative_eta

    subroutine Calculate_Tri_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:) :: D_Matrix

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points
        
        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(dSFMatT(Num_Gauss_Points,Properties%Elements(i)%Number_of_Nodes,2), dSFMat(Num_Gauss_Points,2,Properties%Elements(i)%Number_of_Nodes))

        call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        D_Matrix = 0.0_8

        do k = 1, Properties%Elements(i)%Number_of_Nodes

            do j = 1, Num_Gauss_Points

                dSFMat(j,1,k) = Generate_Tri_Shape_Functions_Derivative_xi(N%Degree, k, eta(j), xi(j))
                dSFMat(j,2,k) = Generate_Tri_Shape_Functions_Derivative_eta(N%Degree, k, eta(j), xi(j))

                dSFMatT(j,k,1) = dSFMat(j,1,k)
                dSFMatT(j,k,2) = dSFMat(j,2,k)

            end do

        end do

        do j = 1, Num_Gauss_Points

            D_Matrix = D_Matrix + 0.5_8*w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), (Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do

        if (Properties%g == 1) D_Matrix = D_Matrix*SUM(Properties%Elements(i)%Coordinates(:,1))/size(Properties%Elements(i)%Coordinates,1)

    end subroutine Calculate_Tri_D_Matrix

    subroutine Calculate_Tri_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, Num_Gauss_Points, index_1, index_2, Num_Nodes

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Mass_Matrix = 0.0_8

        do index_1 = 1, Num_Nodes

            do index_2 = 1, Num_Nodes

                do j = 1, Num_Gauss_Points

                    Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + 0.5_8*w(j)*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))*Generate_Triangular_Shape_Functions(N%Degree, index_2, eta(j), xi(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

                end do

            end do

        end do

        if (Properties%g == 1) Mass_Matrix = Mass_Matrix*SUM(Properties%Elements(i)%Coordinates(:,1))/size(Properties%Elements(i)%Coordinates,1)

    end subroutine Calculate_Tri_Mass_Matrix

    subroutine Calculate_Tri_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, Num_Gauss_Points, index_1, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Source_Vector = 0.0_8

        do index_1 = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                Source_Vector(index_1) = Source_Vector(index_1) + 0.5_8*w(j)*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

            end do

        end do

        if (Properties%g == 1) Source_Vector = Source_Vector*SUM(Properties%Elements(i)%Coordinates(:,1))/size(Properties%Elements(i)%Coordinates,1)

    end subroutine Calculate_Tri_Source_Vector

    Function Tri_Side_Boundary(Properties,N,i,Side_Nodes) Result(B_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer, dimension(:), intent(inout) :: Side_Nodes
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p
        real(kind = 8), dimension(:,:), allocatable :: B_Matrix
        real(kind = 8), dimension(:), allocatable :: xi, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable  :: dx_dxi, dy_dxi

        allocate(B_Matrix(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        Num_Gauss_Points = (N%Degree+1)

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(Nodes(N%Degree+1))

        call Generate_1D_Triangular_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do l = 3, N%Degree+1
            if (l == 3) then
                Nodes(l) = 4
            else
                Nodes(l) = 1 + Nodes(l-1)
            end if
        end do

        if (Side_Nodes(2) == 3 .and. Side_Nodes(1) == 1) then
            Side_Nodes(1) = 3
            Side_Nodes(2) = 1
        end if

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8

        Shape_Functions = 0.0_8

        do l = 1, Num_Gauss_Points

            do p = 1, size(Nodes)
                k = Side_Nodes(p)
                Shape_Functions(k,l) = Generate_Triangular_Shape_Functions(N%Degree, Nodes(p), 0.0_8, xi(l))
                dx_dxi(l) = dx_dxi(l) + Generate_Tri_Shape_Functions_Derivative_xi(N%Degree, Nodes(p), 0.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Tri_Shape_Functions_Derivative_xi(N%Degree, Nodes(p), 0.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
            end do

        end do

        B_Matrix = 0.0_8

        do a = 1, size(Side_Nodes)

           do b = 1, size(Side_Nodes)

                index_1 = Side_Nodes(a)
                index_2 = Side_Nodes(b)

                B_Matrix(index_1,index_2) = B_Matrix(index_1,index_2) + 0.5_8*sum(w*sqrt(dx_dxi**2 + dy_dxi**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                
            end do

        end do

    end Function Tri_Side_Boundary

end module m_Create_Triangular_Shape_Functions_D