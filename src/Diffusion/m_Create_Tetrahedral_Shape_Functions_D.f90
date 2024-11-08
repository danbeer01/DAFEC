module m_Create_Tetrahedral_Shape_Functions_D
!
! Purpose:
! To create the shape functions for a Tetrahedral element
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

    Function Generate_Tetrahedral_Shape_Functions(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: eta, xi, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions(4))

            Shape_Functions(1) = 1.0_8 - xi - eta - zeta
            Shape_Functions(2) = xi
            Shape_Functions(3) = eta
            Shape_Functions(4) = zeta

            Value = Shape_Functions(a)

        else if (Degree == 2) then

            allocate(Shape_Functions(10))

            Shape_Functions(1) = (1.0_8 - xi - eta - zeta)*(2.0_8*(1.0_8 - xi - eta - zeta) - 1.0_8)
            Shape_Functions(2) = xi*(2.0_8*xi - 1.0_8)
            Shape_Functions(3) = eta*(2.0_8*eta - 1.0_8)
            Shape_Functions(4) = zeta*(2.0_8*zeta - 1.0_8)
            Shape_Functions(5) = 4.0_8*xi*(1.0_8 - xi - eta - zeta)
            Shape_Functions(6) = 4.0_8*xi*eta
            Shape_Functions(7) = 4.0_8*eta*(1.0_8 - xi - eta - zeta)
            Shape_Functions(8) = 4.0_8*zeta*(1.0_8 - xi - eta - zeta)
            Shape_Functions(9) = 4.0_8*zeta*xi
            Shape_Functions(10) = 4.0_8*zeta*eta

            Value = Shape_Functions(a)

        end if

    end Function Generate_Tetrahedral_Shape_Functions

    Function Generate_Tetrahedral_Shape_Functions_Derivative_xi(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(4))

            Shape_Functions_Derivative(1) = -1.0_8
            Shape_Functions_Derivative(2) = 1.0_8
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 0.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(10))

            Shape_Functions_Derivative(1) = -(2.0_8*(1.0_8 - xi - eta - zeta) - 1.0_8) - 2.0_8*(1.0_8 - xi - eta - zeta)
            Shape_Functions_Derivative(2) = 2.0_8*xi - 1.0_8 + 2.0_8*xi
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 0.0_8
            Shape_Functions_Derivative(5) = 4.0_8*(1.0_8 - xi - eta - zeta) - 4.0_8*xi
            Shape_Functions_Derivative(6) = 4.0_8*eta
            Shape_Functions_Derivative(7) = -4.0_8*eta
            Shape_Functions_Derivative(8) = -4.0_8*zeta
            Shape_Functions_Derivative(9) = 4.0_8*zeta
            Shape_Functions_Derivative(10) = 0.0_8

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Tetrahedral_Shape_Functions_Derivative_xi

    Function Generate_Tetrahedral_Shape_Functions_Derivative_eta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(4))

            Shape_Functions_Derivative(1) = -1.0_8
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 1.0_8
            Shape_Functions_Derivative(4) = 0.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(10))

            Shape_Functions_Derivative(1) = -(2.0_8*(1.0_8 - xi - eta - zeta) - 1.0_8) - 2.0_8*(1.0_8 - xi - eta - zeta)
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 2.0_8*eta - 1.0_8 + 2.0_8*eta
            Shape_Functions_Derivative(4) = 0.0_8
            Shape_Functions_Derivative(5) = -4.0_8*xi
            Shape_Functions_Derivative(6) = 4.0_8*xi
            Shape_Functions_Derivative(7) = 4.0_8*(1.0_8 - xi - eta - zeta) - 4.0_8*eta
            Shape_Functions_Derivative(8) = -4.0_8*zeta
            Shape_Functions_Derivative(9) = 0.0_8
            Shape_Functions_Derivative(10) = 4.0_8*zeta

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Tetrahedral_Shape_Functions_Derivative_eta

    Function Generate_Tetrahedral_Shape_Functions_Derivative_zeta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(4))

            Shape_Functions_Derivative(1) = -1.0_8
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 1.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(10))

            Shape_Functions_Derivative(1) = -(2.0_8*(1.0_8 - xi - eta - zeta) - 1.0_8) - 2.0_8*(1.0_8 - xi - eta - zeta)
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(4) = 2.0_8*zeta - 1.0_8 + 2.0_8*zeta
            Shape_Functions_Derivative(5) = -4.0_8*xi
            Shape_Functions_Derivative(6) = 0.0_8
            Shape_Functions_Derivative(7) = -4.0_8*eta
            Shape_Functions_Derivative(8) = 4.0_8*(1.0_8 - xi - eta - zeta) - 4.0_8*zeta
            Shape_Functions_Derivative(9) = 4.0_8*xi
            Shape_Functions_Derivative(10) = 4.0_8*eta

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Tetrahedral_Shape_Functions_Derivative_zeta

    subroutine Calculate_Tet_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:)   :: D_Matrix

        if (N%Degree == 1) then
            Num_Gauss_Points = 5
        else if (N%Degree == 2) then
            Num_Gauss_Points = 11
        end if

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(dSFMatT(Num_Gauss_Points,Num_Nodes,3), dSFMat(Num_Gauss_Points,3,Num_Nodes))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Tetrahedral_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        D_Matrix = 0.0_8

        do k = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                dSFMat(j,1,k) = Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,2,k) = Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,3,k) = Generate_Tetrahedral_Shape_Functions_Derivative_zeta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))

                dSFMatT(j,k,1) = dSFMat(j,1,k)
                dSFMatT(j,k,2) = dSFMat(j,2,k)
                dSFMatT(j,k,3) = dSFMat(j,3,k)

            end do

        end do

        do j = 1, Num_Gauss_Points

            D_Matrix = D_Matrix + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(matmul(dSFMatT(j,:,:), Properties%Elements(i)%Inverse_Jacobian(j,:,:)), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do
        
    end subroutine Calculate_Tet_D_Matrix

    subroutine Calculate_Tet_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, index_2, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        if (N%Degree == 1) then
            Num_Gauss_Points = 5
        else if (N%Degree == 2) then
            Num_Gauss_Points = 11
        end if

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Tetrahedral_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Mass_Matrix = 0.0_8

        do index_1 = 1, Num_Nodes

            do index_2 = 1, Num_Nodes

                do j = 1, Num_Gauss_Points

                    Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*(Generate_Tetrahedral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*Generate_Tetrahedral_Shape_Functions(Properties, N%Degree, index_2, eta(j), xi(j), zeta(j)))*ABS(Properties%Elements(i)%Det_Jacobian(j))

                end do

            end do

        end do

    end subroutine Calculate_Tet_Mass_Matrix

    subroutine Calculate_Tet_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        if (N%Degree == 1) then
            Num_Gauss_Points = 5
        else if (N%Degree == 2) then
            Num_Gauss_Points = 11
        end if

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Tetrahedral_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Source_Vector = 0.0_8

        do index_1 = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                Source_Vector(index_1) = Source_Vector(index_1) + w(j)*(Generate_Tetrahedral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j)))*ABS(Properties%Elements(i)%Det_Jacobian(j))

            end do

        end do

    end subroutine Calculate_Tet_Source_Vector

    Function Tet_Face_Boundary(Properties,N,i,side_nodes) Result(B_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer, dimension(:), intent(inout) :: side_nodes
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p
        real(kind = 8), dimension(:,:), allocatable :: B_Matrix
        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi, dz_dxi, dx_deta, dy_deta, dz_deta

        allocate(B_Matrix(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        if(N%Degree == 1) then
            Num_Gauss_Points = 3
        else
            Num_Gauss_Points = 7
        end if

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(Nodes(((N%Degree+1)*(N%Degree+2)/2)))

        call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        Nodes(3) = 3
        if(N%Degree == 2) then
            Nodes(4:6) = (/5, 6, 7/)
        end if

        if (N%Degree == 2) then
            if (side_nodes(1) == 1 .and. side_nodes(2) == 2 .and. side_nodes(3) == 4) side_nodes = (/1,2,4,5,9,8/)
            if (side_nodes(1) == 2 .and. side_nodes(2) == 3 .and. side_nodes(3) == 4) side_nodes = (/2,3,4,6,10,9/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 3 .and. side_nodes(3) == 4) side_nodes = (/1,3,4,7,10,8/)
        end if

        B_Matrix = 0.0_8

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        dz_dxi = 0.0_8
        dx_deta = 0.0_8
        dy_deta = 0.0_8
        dz_deta = 0.0_8

        Shape_Functions = 0.0_8

        do l = 1, Num_Gauss_Points

            do p = 1, size(Nodes)
                k = side_nodes(p)
                Shape_Functions(k,l) = Generate_Tetrahedral_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)
                dx_dxi(l) = dx_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                dz_dxi(l) = dz_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
                dx_deta(l) = dx_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                dy_deta(l) = dy_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                dz_deta(l) = dz_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
            end do

        end do

        do a = 1, size(side_nodes)

           do b = 1, size(side_nodes)

                index_1 = side_nodes(a)
                index_2 = side_nodes(b)

                B_Matrix(index_1,index_2) = B_Matrix(index_1,index_2) + 0.25_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*(Shape_Functions(index_1,:)*Shape_Functions(index_2,:)))

            end do

        end do

    end Function Tet_Face_Boundary

end module