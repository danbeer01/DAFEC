module m_Create_Prismatic_Shape_Functions
!
! Purpose:
! To create the shape functions for a Prismatic element
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties
use m_Gauss_Points

implicit none

contains

    Function Generate_Prismatic_Shape_Functions(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, t
        real(kind = 8), intent(in)                :: eta, xi, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions(6))

            Shape_Functions(1) = (1.0_8 - xi - eta)*(1.0_8 - zeta)
            Shape_Functions(2) = xi*(1.0_8 - zeta)
            Shape_Functions(3) = eta*(1.0_8 - zeta)
            Shape_Functions(4) = (1.0_8 - xi - eta)*(1.0_8 + zeta)
            Shape_Functions(5) = xi*(1.0_8 + zeta)
            Shape_Functions(6) = eta*(1.0_8 + zeta)

            Value = Shape_Functions(a)*0.5_8

        else if (Degree == 2) then

            allocate(Shape_Functions(15))

            t = 1.0_8 - xi - eta

            Shape_Functions(1) = t*((2.0_8*t - 1.0_8)*(1.0_8 - zeta) - (1.0_8 - zeta**2))
            Shape_Functions(2) = xi*((2.0_8*xi - 1.0_8)*(1.0_8 - zeta)  - (1.0_8 - zeta**2))
            Shape_Functions(3) = eta*((2.0_8*eta - 1.0_8)*(1.0_8 - zeta)  - (1.0_8 - zeta**2))

            Shape_Functions(4) = t*((2.0_8*t - 1.0_8)*(1.0_8 + zeta)  - (1.0_8 - zeta**2))
            Shape_Functions(5) = xi*((2.0_8*xi - 1.0_8)*(1.0_8 + zeta)  - (1.0_8 - zeta**2))
            Shape_Functions(6) = eta*((2.0_8*eta - 1.0_8)*(1.0_8 + zeta)  - (1.0_8 - zeta**2))

            Shape_Functions(7) = 4.0_8*t*xi*(1.0_8 - zeta)
            Shape_Functions(8) = 4.0_8*xi*eta*(1.0_8 - zeta)
            Shape_Functions(9) = 4.0_8*eta*t*(1.0_8 - zeta)

            Shape_Functions(10) = 2.0_8*t*(1.0_8 - zeta**2)
            Shape_Functions(11) = 2.0_8*xi*(1.0_8 - zeta**2)
            Shape_Functions(12) = 2.0_8*eta*(1.0_8 - zeta**2)

            Shape_Functions(13) = 4.0_8*t*xi*(1.0_8 + zeta)
            Shape_Functions(14) = 4.0_8*xi*eta*(1.0_8 + zeta)
            Shape_Functions(15) = 4.0_8*eta*t*(1.0_8 + zeta)

            Value = Shape_Functions(a)*0.5_8

        end if

    end Function Generate_Prismatic_Shape_Functions

    Function Generate_Prismatic_Shape_Functions_Derivative_xi(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, t
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(6))

            Shape_Functions_Derivative(1) = -1.0_8 + zeta
            Shape_Functions_Derivative(2) = 1.0_8 - zeta
            Shape_Functions_Derivative(3) = 0.0_8
            Shape_Functions_Derivative(5) = (1.0_8 + zeta)
            Shape_Functions_Derivative(6) = 0.0_8
            Shape_Functions_Derivative(4) = -(1.0_8 + zeta)

            Value = Shape_Functions_Derivative(a)*0.5_8

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(15))

            t = 1.0_8 - xi - eta

            Shape_Functions_Derivative(1) = -(2.0_8*t - 1.0_8)*(1.0_8 - zeta) + (1.0_8 - zeta**2) - 2.0_8*t*(1.0_8 - zeta)
            Shape_Functions_Derivative(2) = (2.0_8*xi - 1.0_8)*(1.0_8 - zeta)  - (1.0_8 - zeta**2) + 2.0_8*xi*(1.0_8 - zeta)
            Shape_Functions_Derivative(3) = 0.0_8

            Shape_Functions_Derivative(4) = -(2.0_8*t - 1.0_8)*(1.0_8 + zeta)  + (1.0_8 - zeta**2) - 2.0_8*t*(1.0_8 + zeta)
            Shape_Functions_Derivative(5) = (2.0_8*xi - 1.0_8)*(1.0_8 + zeta)  - (1.0_8 - zeta**2) + 2.0_8*xi*(1.0_8 + zeta)
            Shape_Functions_Derivative(6) = 0.0_8

            Shape_Functions_Derivative(7) = 4.0_8*t*(1.0_8 - zeta) - 4.0_8*xi*(1.0_8 - zeta)
            Shape_Functions_Derivative(8) = 4.0_8*(1.0_8 - zeta)*eta
            Shape_Functions_Derivative(9) = -4.0_8*(1.0_8 - zeta)*eta

            Shape_Functions_Derivative(10) = -2.0_8*(1.0_8 - zeta**2)
            Shape_Functions_Derivative(11) = 2.0_8*(1.0_8 - zeta**2)
            Shape_Functions_Derivative(12) = 0.0_8

            Shape_Functions_Derivative(13) = 4.0_8*t*(1.0_8 + zeta) - 4.0_8*xi*(1.0_8 + zeta)
            Shape_Functions_Derivative(14) = 4.0_8*(1.0_8 + zeta)*eta
            Shape_Functions_Derivative(15) = -4.0_8*(1.0_8 + zeta)*eta

            Value = Shape_Functions_Derivative(a)*0.5_8

        end if

    end Function Generate_Prismatic_Shape_Functions_Derivative_xi

    Function Generate_Prismatic_Shape_Functions_Derivative_eta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, t
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(6))

            Shape_Functions_Derivative(1) = -1.0_8 + zeta
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = 1.0_8 - zeta
            Shape_Functions_Derivative(4) = -(1.0_8 + zeta)
            Shape_Functions_Derivative(5) = 0.0_8
            Shape_Functions_Derivative(6) = (1.0_8 + zeta)

            Value = Shape_Functions_Derivative(a)*0.5_8

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(15))

            t = 1.0_8 - xi - eta

            Shape_Functions_Derivative(1) = -(2.0_8*t - 1.0_8)*(1.0_8 - zeta) + (1.0_8 - zeta**2) - 2.0_8*t*(1.0_8 - zeta)
            Shape_Functions_Derivative(2) = 0.0_8
            Shape_Functions_Derivative(3) = (2.0_8*eta - 1.0_8)*(1.0_8 - zeta)  - (1.0_8 - zeta**2) + 2.0_8*eta*(1.0_8 - zeta)

            Shape_Functions_Derivative(4) = -(2.0_8*t - 1.0_8)*(1.0_8 + zeta)  + (1.0_8 - zeta**2) - 2.0_8*t*(1.0_8 + zeta)
            Shape_Functions_Derivative(5) = 0.0_8
            Shape_Functions_Derivative(6) = (2.0_8*eta - 1.0_8)*(1.0_8 + zeta)  - (1.0_8 - zeta**2) + 2.0_8*eta*(1.0_8 + zeta)

            Shape_Functions_Derivative(7) = -4.0_8*xi*(1.0_8 - zeta)
            Shape_Functions_Derivative(8) = 4.0_8*(1.0_8 - zeta)*xi
            Shape_Functions_Derivative(9) = -4.0_8*(1.0_8 - zeta)*eta + 4.0_8*t*(1.0_8 - zeta)

            Shape_Functions_Derivative(10) = -2.0_8*(1.0_8 - zeta**2)
            Shape_Functions_Derivative(11) = 0.0_8
            Shape_Functions_Derivative(12) = 2.0_8*(1.0_8 - zeta**2)

            Shape_Functions_Derivative(13) = -4.0_8*xi*(1.0_8 + zeta)
            Shape_Functions_Derivative(14) = 4.0_8*(1.0_8 + zeta)*xi
            Shape_Functions_Derivative(15) = 4.0_8*t*(1.0_8 + zeta) -4.0_8*(1.0_8 + zeta)*eta

            Value = Shape_Functions_Derivative(a)*0.5_8

        end if

    end Function Generate_Prismatic_Shape_Functions_Derivative_eta

    Function Generate_Prismatic_Shape_Functions_Derivative_zeta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, t
        real(kind = 8), intent(in)                :: xi, eta, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(6))

            Shape_Functions_Derivative(1) = -1.0_8 + xi + eta
            Shape_Functions_Derivative(2) = -xi
            Shape_Functions_Derivative(3) = -eta
            Shape_Functions_Derivative(4) = 1.0_8 - xi - eta
            Shape_Functions_Derivative(5) = xi
            Shape_Functions_Derivative(6) = eta

            Value = Shape_Functions_Derivative(a)*0.5_8

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(15))

            t = 1.0_8 - xi - eta

            Shape_Functions_Derivative(1) = -t*(2.0_8*t - 1.0_8) + 2.0_8*zeta*t
            Shape_Functions_Derivative(2) = -xi*(2.0_8*xi - 1.0_8) + 2.0_8*zeta*xi
            Shape_Functions_Derivative(3) = -eta*(2.0_8*eta - 1.0_8) + 2.0_8*zeta*eta

            Shape_Functions_Derivative(4) = t*(2.0_8*t - 1.0_8) + 2.0_8*zeta*t
            Shape_Functions_Derivative(5) = xi*(2.0_8*xi - 1.0_8) + 2.0_8*zeta*xi
            Shape_Functions_Derivative(6) = eta*(2.0_8*eta - 1.0_8) + 2.0_8*zeta*eta

            Shape_Functions_Derivative(7) = -4.0_8*xi*t
            Shape_Functions_Derivative(8) = -4.0_8*xi*eta
            Shape_Functions_Derivative(9) = -4.0_8*eta*t

            Shape_Functions_Derivative(10) = -4.0_8*zeta*t
            Shape_Functions_Derivative(11) = -4.0_8*zeta*xi
            Shape_Functions_Derivative(12) = -4.0_8*zeta*eta

            Shape_Functions_Derivative(13) = 4.0_8*xi*t
            Shape_Functions_Derivative(14) = 4.0_8*xi*eta
            Shape_Functions_Derivative(15) = 4.0_8*eta*t

            Value = Shape_Functions_Derivative(a)*0.5_8

        end if

    end Function Generate_Prismatic_Shape_Functions_Derivative_zeta

    subroutine Calculate_Pris_Streaming_Matrix(Properties,N,i,Streaming_Matrix,mu_ang,eta_ang,xi_ang)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points

        real(kind = 8), intent(in)                :: mu_ang, eta_ang, xi_ang

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:) :: Streaming_Matrix

        if (N%Degree == 1) then

            Num_Gauss_Points = 8

        else if (N%Degree == 2) then

            Num_Gauss_Points = 21

        end if

        allocate(dSFMatT(Num_Gauss_Points,Properties%Elements(i)%Number_of_Nodes,3), dSFMat(Num_Gauss_Points,3,Properties%Elements(i)%Number_of_Nodes))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Prismatic_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Streaming_Matrix = 0.0_8

        do j = 1, Num_Gauss_Points

            do k = 1, Properties%Elements(i)%Number_of_Nodes

                dSFMat(j,1,k) = mu_ang*Generate_Prismatic_Shape_Functions(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,2,k) = eta_ang*Generate_Prismatic_Shape_Functions(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,3,k) = xi_ang*Generate_Prismatic_Shape_Functions(Properties, N%Degree, k, eta(j), xi(j), zeta(j))

                dSFMatT(j,k,1) = Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMatT(j,k,2) = Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMatT(j,k,3) = Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))

            end do

            Streaming_Matrix = Streaming_Matrix + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(dSFMatT(j,:,:), (Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do
        
    end subroutine Calculate_Pris_Streaming_Matrix

    subroutine Calculate_Pris_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: index_1, index_2, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        if (N%Degree == 1) then

            Num_Gauss_Points = 8

        else if (N%Degree == 2) then

            Num_Gauss_Points = 21

        end if

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Prismatic_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Mass_Matrix = 0.0_8

        do index_1 = 1, Num_Nodes

            do index_2 = 1, Num_Nodes

                do j = 1, Num_Gauss_Points

                    Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Prismatic_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*Generate_Prismatic_Shape_Functions(Properties, N%Degree, index_2, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

                end do

            end do

        end do

    end subroutine Calculate_Pris_Mass_Matrix

    subroutine Calculate_Pris_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: index_1, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        if (N%Degree == 1) then

            Num_Gauss_Points = 8

        else if (N%Degree == 2) then

            Num_Gauss_Points = 21

        end if
        
        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Prismatic_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Source_Vector = 0.0_8

        do index_1 = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Prismatic_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

            end do

        end do

    end subroutine Calculate_Pris_Source_Vector

    Function Integrate_Pris_Face(Properties,N,i,j) Result(F_out)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p
        real(kind = 8), dimension(:,:), allocatable :: F_out
        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi, dz_dxi, dx_deta, dy_deta, dz_deta

        allocate(F_out(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        if (j == 1 .or. j == 2 .or. j == 3) then

            Num_Gauss_Points = (N%Degree+1)**2

            allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))
    
            allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))
    
            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))
    
            allocate(Nodes((N%Degree)*4))
    
            call Generate_2D_Pris_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            Nodes(1) = 1
            Nodes(2) = 2
            Nodes(3) = 5
            Nodes(4) = 4
            if (N%Degree == 2) then
                Nodes(5) = 7
                Nodes(6) = 11
                Nodes(7) = 13
                Nodes(8) = 10
            end if

            dx_dxi = 0.0_8
            dy_dxi = 0.0_8
            dz_dxi = 0.0_8
            dx_deta = 0.0_8
            dy_deta = 0.0_8
            dz_deta = 0.0_8

            Shape_Functions = 0.0_8

            do l = 1, Num_Gauss_Points

                do p = 1, size(Nodes)
                    k = Properties%Elements(i)%Side_Nodes(j,p)
                    Shape_Functions(k,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))
                    dx_dxi(l) = dx_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,3)
                end do

            end do

        else

            if(N%Degree == 1) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if

            allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

            allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))

            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

            allocate(Nodes(((N%Degree+1)*(N%Degree+2)/2)))

            call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            Nodes(1) = 1
            Nodes(2) = 2
            Nodes(3) = 3
            if (N%Degree == 2) then
                Nodes(4) = 7
                Nodes(5) = 8
                Nodes(6) = 9
            end if

            dx_dxi = 0.0_8
            dy_dxi = 0.0_8
            dz_dxi = 0.0_8
            dx_deta = 0.0_8
            dy_deta = 0.0_8
            dz_deta = 0.0_8
    
            Shape_Functions = 0.0_8
  
            do l = 1, Num_Gauss_Points
    
                do p = 1, size(Nodes)
                    k = Properties%Elements(i)%Side_Nodes(j,p)
                    Shape_Functions(k,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)
                    dx_dxi(l) = dx_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,3)
                end do
    
            end do

        end if

        F_out = 0.0_8

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

           do b = 1,  size(Properties%Elements(i)%Side_Nodes(j,:))

                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(i)%Side_Nodes(j,b)

                if (j==1 .or. j==2 .or. j==3) then
                    F_out(index_1,index_2) = F_out(index_1,index_2) + sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                else
                    if (index_1 /= 0 .and. index_2 /= 0) then
                        F_out(index_1,index_2) = F_out(index_1,index_2) + 0.5_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                    end if
                end if

            end do

        end do

    end Function Integrate_Pris_Face

    Function Integrate_Pris_Face_F_in(Properties,N,i,j) Result(F_in)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j
        integer             :: Num_Gauss_Points, a, b, k, k_n, l, p, index_1, index_2

        real(kind = 8), dimension(:,:), allocatable :: F_in
        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        integer, dimension(:), allocatable :: Nodes, Neighbour_Nodes

        real(kind = 8), dimension(:), allocatable :: dx_dxi, dy_dxi, dz_dxi, dx_deta, dy_deta, dz_deta

        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions, Shape_Functions_Neighbour

        allocate(F_in(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes))

        if (j==1 .or. j==2 .or. j==3) then

            Num_Gauss_Points = (N%Degree+1)**2

            allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

            allocate(Shape_Functions_Neighbour(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes,Num_Gauss_Points))
    
            allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))
    
            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))
    
            allocate(Nodes((N%Degree)*4))

            allocate(Neighbour_Nodes((N%Degree)*4))
    
            call Generate_2D_Pris_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            Nodes(1) = 1
            Nodes(2) = 2
            Nodes(3) = 5
            Nodes(4) = 4
            if (N%Degree == 2) then
                Nodes(5) = 7
                Nodes(6) = 11
                Nodes(7) = 13
                Nodes(8) = 10
            end if            

            Neighbour_Nodes = Properties%Elements(i)%Sides(j)%Neighbour_Nodes

            dx_dxi = 0.0_8
            dy_dxi = 0.0_8
            dz_dxi = 0.0_8
            dx_deta = 0.0_8
            dy_deta = 0.0_8
            dz_deta = 0.0_8
    
            Shape_Functions = 0.0_8
    
            Shape_Functions_Neighbour = 0.0_8
    
            do l = 1, Num_Gauss_Points
    
                do p = 1, size(Nodes)
                    k = Properties%Elements(i)%Side_Nodes(j,p)
                    Shape_Functions(k,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))
                    k_n = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),p)
                    Shape_Functions_Neighbour(k_n,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Neighbour_Nodes(p), 0.0_8, xi(l), eta(l))
                    dx_dxi(l) = dx_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_zeta(Properties, N%Degree, Nodes(p), 0.0_8, xi(l), eta(l))*Properties%Elements(i)%Coordinates(k,3)
                end do
    
            end do

        else

            if(N%Degree == 1) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if
    
            allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

            allocate(Shape_Functions_Neighbour(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes,Num_Gauss_Points))
    
            allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))
    
            allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

            allocate(Nodes(((N%Degree+1)*(N%Degree+2)/2)))

            allocate(Neighbour_Nodes(((N%Degree+1)*(N%Degree+2)/2)))

            call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

            Nodes(1) = 1
            Nodes(2) = 2
            Nodes(3) = 3
            if (N%Degree == 2) then
                Nodes(4) = 7
                Nodes(5) = 8
                Nodes(6) = 9
            end if

            Neighbour_Nodes = Properties%Elements(i)%Sides(j)%Neighbour_Nodes

            dx_dxi = 0.0_8
            dy_dxi = 0.0_8
            dz_dxi = 0.0_8
            dx_deta = 0.0_8
            dy_deta = 0.0_8
            dz_deta = 0.0_8
    
            Shape_Functions = 0.0_8
    
            Shape_Functions_Neighbour = 0.0_8
  
            do l = 1, Num_Gauss_Points
    
                do p = 1, size(Nodes)
                    k = Properties%Elements(i)%Side_Nodes(j,p)
                    Shape_Functions(k,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)
                    k_n = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),p)
                    Shape_Functions_Neighbour(k_n,l) = Generate_Prismatic_Shape_Functions(Properties, N%Degree, Neighbour_Nodes(p), eta(l), xi(l), -1.0_8)
                    dx_dxi(l) = dx_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Prismatic_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Prismatic_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), -1.0_8)*Properties%Elements(i)%Coordinates(k,3)
                end do
    
            end do

        end if

        F_in = 0.0_8

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

            do b = 1, size(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(j,:))
 
                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),b)

                if (j==1 .or. j==2 .or. j==3) then
                    F_in(index_1,index_2) = F_in(index_1,index_2) + sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions_Neighbour(index_2,:))
                else
                    if (index_1 /= 0 .and. index_2 /= 0) then
                        F_in(index_1,index_2) = F_in(index_1,index_2) + 0.5_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions_Neighbour(index_2,:))
                    end if
                end if

            end do

        end do

    end Function Integrate_Pris_Face_F_in

end module m_Create_Prismatic_Shape_Functions