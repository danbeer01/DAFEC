module m_Create_Quadrilateral_Shape_Functions
!
! Purpose:
! To create the shape functions for a quadrilateral element
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

    subroutine Calculate_Isoparametric_Quadrilateral_Coordinates(Degree, this)

        type(PropertiesType), intent(inout) :: this

        real(kind = 8)                              :: x, y
        integer, intent(in)                         :: Degree
        integer :: i, j, k, Coord_index

        allocate(this%Isoparametric_Coordinates((Degree+1)*(Degree+1),2))

        Coord_index = 1

        k = 0

        do j = 1,Degree,2

            do i = 1,4

                x = 1.0_8 - 2.0_8*k/real(Degree,kind=8)
                y = 1.0_8 - 2.0_8*k/real(Degree,kind=8)

                if (i == 1) then
                    x = ABS(x)
                    y = ABS(y)
                else if (i == 2) then
                    x = -ABS(x)
                    y = ABS(y)
                else if (i == 3) then
                    x = -ABS(x)
                    y = -ABS(y)
                else if (i == 4) then
                    x = ABS(x)
                    y = -ABS(y)
                end if
                
                this%Isoparametric_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

            do i = 1,Degree-1-2*k

                x = (1.0_8 - 2.0_8*k/real(Degree,kind=8)) - 2.0_8*i*(1.0_8 - 2.0_8*k/real(Degree,kind=8))/(real(Degree,kind=8)-real(2*k,kind=8))
                y = 1.0_8 - 2.0_8*k/real(Degree,kind=8)

                this%Isoparametric_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

            do i = 1,Degree-1-2*k

                y = (1.0_8 - 2.0_8*k/real(Degree,kind=8)) - 2.0_8*i*(1.0_8 - 2.0_8*k/real(Degree,kind=8))/(real(Degree,kind=8)-real(2*k,kind=8))
                x = -1.0_8 + 2.0_8*k/real(Degree,kind=8)

                this%Isoparametric_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

            do i = 1,Degree-1-2*k

                x = (-1.0_8 + 2.0_8*k/real(Degree,kind=8)) + 2.0_8*i*(1.0_8 - 2.0_8*k/real(Degree,kind=8))/(real(Degree,kind=8)-real(2*k,kind=8))
                y = -1.0_8 + 2.0_8*k/real(Degree,kind=8)

                this%Isoparametric_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

            do i = 1,Degree-1-2*k

                y = (-1.0_8 + 2.0_8*k/real(Degree,kind=8)) + 2.0_8*i*(1.0_8 - 2.0_8*k/real(Degree,kind=8))/(real(Degree,kind=8)-real(2*k,kind=8))
                x = 1.0_8 - 2.0_8*k/real(Degree,kind=8)

                this%Isoparametric_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

            k = k + 1

        end do

        if (mod(Degree,2) == 0) this%Isoparametric_Coordinates(Coord_index,:) = [0.0_8, 0.0_8]

    end subroutine Calculate_Isoparametric_Quadrilateral_Coordinates

    Function Generate_Quadrilateral_Shape_Functions(this, Degree, a, eta, xi) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: eta, xi

        integer :: i

        allocate(Shape_Functions((Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions(1) = 0.25_8*(1.0_8 + xi)*(1.0_8 + eta)
            Shape_Functions(2) = 0.25_8*(1.0_8 - xi)*(1.0_8 + eta)
            Shape_Functions(3) = 0.25_8*(1.0_8 - xi)*(1.0_8 - eta)
            Shape_Functions(4) = 0.25_8*(1.0_8 + xi)*(1.0_8 - eta)

            Value = Shape_Functions(a)

        else if (Degree == 2) then

            Shape_Functions(1) = 0.25_8*xi*(1.0_8 + xi)*eta*(1.0_8 + eta)
            Shape_Functions(2) = -0.25_8*xi*(1.0_8 - xi)*eta*(1.0_8 + eta)
            Shape_Functions(3) = 0.25_8*xi*(1.0_8 - xi)*eta*(1.0_8 - eta)
            Shape_Functions(4) = -0.25_8*xi*(1.0_8 + xi)*eta*(1.0_8 - eta)
            Shape_Functions(5) = 0.5_8*(1.0_8 - xi*xi)*eta*(1.0_8 + eta)
            Shape_Functions(6) = -0.5_8*xi*(1.0_8 - xi)*(1.0_8 - eta*eta)
            Shape_Functions(7) = -0.5_8*(1.0_8 - xi*xi)*eta*(1.0_8 - eta)
            Shape_Functions(8) = 0.5_8*xi*(1.0_8 + xi)*(1.0_8 - eta*eta)
            Shape_Functions(9) = (1.0_8 - xi*xi)*(1.0_8 - eta*eta)

            Value = Shape_Functions(a)

        else

            Value = 1.0_8

            do i = 1,(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        Value = Value*(xi - this%Isoparametric_Coordinates(i,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6) then

                        Value = Value*(eta - this%Isoparametric_Coordinates(i,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    end if

                end if

            end do

        end if

    end Function Generate_Quadrilateral_Shape_Functions

    Function Generate_Shape_Functions_Derivative_xi(this, Degree, a, eta, xi) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, temp
        real(kind = 8), intent(in)                :: xi, eta

        integer :: i, j

        allocate(Shape_Functions_Derivative((Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 0.25_8*(1.0_8 + eta)
            Shape_Functions_Derivative(2) = -0.25_8*(1.0_8 + eta)
            Shape_Functions_Derivative(3) = -0.25_8*(1.0_8 - eta)
            Shape_Functions_Derivative(4) = 0.25_8*(1.0_8 - eta)

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then
    
            Shape_Functions_Derivative(1) = 0.25_8*eta*(1.0_8 + eta)*(1.0_8 + 2.0_8*xi)
            Shape_Functions_Derivative(2) = -0.25_8*eta*(1.0_8 + eta)*(1.0_8 - 2.0_8*xi)
            Shape_Functions_Derivative(3) = 0.25_8*eta*(1.0_8 - eta)*(1.0_8 - 2.0_8*xi)
            Shape_Functions_Derivative(4) = -0.25_8*eta*(1.0_8 - eta)*(1.0_8 + 2.0_8*xi)
            Shape_Functions_Derivative(5) = 0.5_8*eta*(1.0_8 + eta)*(- 2.0_8*xi)
            Shape_Functions_Derivative(6) = -0.5_8*(1.0_8 - eta*eta)*(1.0_8 - 2.0_8*xi)
            Shape_Functions_Derivative(7) = -0.5_8*eta*(1.0_8 - eta)*(- 2.0_8*xi)
            Shape_Functions_Derivative(8) = 0.5_8*(1.0_8 - eta*eta)*(1.0_8 + 2.0_8*xi)
            Shape_Functions_Derivative(9) = (1.0_8 - eta*eta)*(- 2.0_8*xi)

            Value = Shape_Functions_Derivative(a)

        else

            Value = 0.0_8

            do i = 1,(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) <1e-6) then

                        temp = 1.0_8

                        do j = 1,(Degree+1)*(Degree+1)

                            if (j /= a) then

                                if (ABS(this%Isoparametric_Coordinates(j,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                                    if (j /= i) then

                                        temp = temp*(xi - this%Isoparametric_Coordinates(j,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(j,1))

                                    end if

                                end if

                            end if

                        end do

                        Value = Value + temp/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    end if

                end if

            end do

            do i = 1,(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6) then

                        Value = Value*(eta - this%Isoparametric_Coordinates(i,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    end if

                end if

            end do

        end if

    end Function Generate_Shape_Functions_Derivative_xi

    Function Generate_Shape_Functions_Derivative_eta(this, Degree, a, eta, xi) Result(Value)

        type(PropertiesType), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, temp
        real(kind = 8), intent(in)                :: xi, eta

        integer :: i, j

        allocate(Shape_Functions_Derivative((Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 0.25_8*(1.0_8 + xi)
            Shape_Functions_Derivative(2) = 0.25_8*(1.0_8 - xi)
            Shape_Functions_Derivative(3) = -0.25_8*(1.0_8 - xi)
            Shape_Functions_Derivative(4) = -0.25_8*(1.0_8 + xi)

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            Shape_Functions_Derivative(1) = 0.25_8*xi*(1.0_8 + xi)*(1.0_8 + 2.0_8*eta)
            Shape_Functions_Derivative(2) = -0.25_8*xi*(1.0_8 - xi)*(1.0_8 + 2.0_8*eta)
            Shape_Functions_Derivative(3) = 0.25_8*xi*(1.0_8 - xi)*(1.0_8 - 2.0_8*eta)
            Shape_Functions_Derivative(4) = -0.25_8*xi*(1.0_8 + xi)*(1.0_8 - 2.0_8*eta)
            Shape_Functions_Derivative(5) = 0.5_8*(1.0_8 - xi*xi)*(1.0_8 + 2.0_8*eta)
            Shape_Functions_Derivative(6) = -0.5_8*xi*(1.0_8 - xi)*(- 2.0_8*eta)
            Shape_Functions_Derivative(7) = -0.5_8*(1.0_8 - xi*xi)*(1.0_8 - 2.0_8*eta)
            Shape_Functions_Derivative(8) = 0.5_8*xi*(1.0_8 + xi)*(- 2.0_8*eta)
            Shape_Functions_Derivative(9) = (1.0_8 - xi*xi)*(- 2.0_8*eta)

            Value = Shape_Functions_Derivative(a)

        else

            Value = 0.0_8

            do i = 1,(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6) then

                        temp = 1.0_8

                        do j = 1,(Degree+1)*(Degree+1)

                            if (j /= a) then

                                if (ABS(this%Isoparametric_Coordinates(j,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6) then

                                    if (j /= i) then

                                        temp = temp*(eta - this%Isoparametric_Coordinates(j,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(j,2))

                                    end if

                                end if

                            end if

                        end do

                        Value = Value + temp/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    end if

                end if

            end do

            do i = 1,(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        Value = Value*(xi - this%Isoparametric_Coordinates(i,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    end if

                end if

            end do

        end if

    end Function Generate_Shape_Functions_Derivative_eta

    subroutine Calculate_Quad_Streaming_Matrix(Properties,N,i,Streaming_Matrix,mu_ang,eta_ang)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points

        real(kind = 8), intent(in)                :: mu_ang, eta_ang

        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        real(kind = 8)                            :: xi_val, eta_val, r

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:) :: Streaming_Matrix

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

        allocate(dSFMatT(Num_Gauss_Points,Properties%Elements(i)%Number_of_Nodes,2), dSFMat(Num_Gauss_Points,2,Properties%Elements(i)%Number_of_Nodes))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Streaming_Matrix = 0.0_8

        do j = 1, Num_Gauss_Points

            xi_val = xi(j) 
            eta_val = eta(j)

            do k = 1, Properties%Elements(i)%Number_of_Nodes

                dSFMat(j,1,k) = mu_ang*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, k, eta(j), xi(j))
                dSFMat(j,2,k) = eta_ang*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, k, eta(j), xi(j))

                dSFMatT(j,k,1) = Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, k, eta(j), xi(j))
                dSFMatT(j,k,2) = Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, k, eta(j), xi(j))

            end do

            r = 0.0_8

            do k = 1, Properties%Elements(i)%Number_of_Nodes

                r = r + Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, k, eta_val, xi_val)*Properties%Elements(i)%Coordinates(k,1)

            end do

            if(Properties%g == 0) Streaming_Matrix = Streaming_Matrix + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))
            if(Properties%g == 1) Streaming_Matrix = Streaming_Matrix + r*w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do

    end subroutine Calculate_Quad_Streaming_Matrix

    subroutine Calculate_Quad_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: index_1, index_2, j, k, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        real(kind = 8)                            :: xi_val, eta_val, r

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Mass_Matrix = 0.0_8

        do j = 1, Num_Gauss_Points

            xi_val = xi(j)

            eta_val = eta(j)

            r = 0.0_8

            do k = 1, Num_Nodes

                r = r + Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, k, eta_val, xi_val)*Properties%Elements(i)%Coordinates(k,1)

            end do

            do index_1 = 1, Num_Nodes

                do index_2 = 1, Num_Nodes

                    if(Properties%g == 0) Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta_val, xi_val)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_2, eta_val, xi_val)*Properties%Elements(i)%Det_Jacobian(j)
                    if(Properties%g == 1) Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + r*w(j)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta_val, xi_val)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_2, eta_val, xi_val)*Properties%Elements(i)%Det_Jacobian(j)

                end do

            end do

        end do

    end subroutine Calculate_Quad_Mass_Matrix

    subroutine Calculate_Quad_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: index_1, j, k, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        real(kind = 8)                            :: xi_val, eta_val, r

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Source_Vector = 0.0_8

        do j = 1, Num_Gauss_Points

            xi_val = xi(j) 

            eta_val = eta(j)

            r = 0.0_8

            do k = 1, Num_Nodes

                r = r + Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, k, eta_val, xi_val)*Properties%Elements(i)%Coordinates(k,1)

            end do

            do index_1 = 1, Num_Nodes

                if(Properties%g == 0) Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta_val, xi_val)*Properties%Elements(i)%Det_Jacobian(j)
                if(Properties%g == 1) Source_Vector(index_1) = Source_Vector(index_1) + r*w(j)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta_val, xi_val)*Properties%Elements(i)%Det_Jacobian(j)

            end do

        end do

    end subroutine Calculate_Quad_Source_Vector

    Function Integrate_Quad_Side(Properties,N,i,j) Result(F_out)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p
        real(kind = 8), dimension(:,:), allocatable :: F_out
        real(kind = 8), dimension(:), allocatable :: xi, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi, r

        allocate(F_out(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        Num_Gauss_Points = (N%Degree+1)

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), r(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(Nodes(N%Degree+1))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do l = 3, N%Degree+1
            if (l == 3) then
                Nodes(l) = 5
            else
                Nodes(l) = 1 + Nodes(l-1)
            end if
        end do

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        r = 0.0_8

        Shape_Functions = 0.0_8

        do l = 1, Num_Gauss_Points

            do p = 1, size(Nodes)
                k = Properties%Elements(i)%Side_Nodes(j,p)
                Shape_Functions(k,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))
                dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
                r(l) = r(l) + Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
            end do

        end do

        F_out = 0.0_8

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

           do b = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(i)%Side_Nodes(j,b)

                if(Properties%g == 0) F_out(index_1,index_2) = F_out(index_1,index_2) + sum(w*sqrt(dx_dxi**2 + dy_dxi**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                if(Properties%g == 1) F_out(index_1,index_2) = F_out(index_1,index_2) + sum(w*sqrt(dx_dxi**2 + dy_dxi**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))

            end do

        end do

    end Function Integrate_Quad_Side

    Function Integrate_Quad_Side_F_in(Properties,N,i,j) Result(F_in)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j
        integer             :: Num_Gauss_Points, a, b, k, k_n, l, p, temp, index_1, index_2

        real(kind = 8), dimension(:,:), allocatable :: F_in
        real(kind = 8), dimension(:), allocatable :: xi, w

        integer, dimension(N%Degree+1) :: Nodes, Neighbour_Nodes

        real(kind = 8), dimension(:), allocatable  :: dx_dxi, dy_dxi, r

        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions, Shape_Functions_Neighbour

        allocate(F_in(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes))

        Num_Gauss_Points = (N%Degree+1)

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes, Num_Gauss_Points))

        allocate(Shape_Functions_Neighbour(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes, Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), r(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do l = 3, N%Degree+1
            if (l == 3) then
                Nodes(l) = 5
            else
                Nodes(l) = 1 + Nodes(l-1)
            end if
        end do

        Neighbour_Nodes = Nodes
        temp = Neighbour_Nodes(1)
        Neighbour_Nodes(1) = Neighbour_Nodes(2)
        Neighbour_Nodes(2) = temp
        do l = 3, 2 + (size(Neighbour_Nodes) - 2)/2
            temp = Neighbour_Nodes(l)
            Neighbour_Nodes(l) = Neighbour_Nodes(size(Neighbour_Nodes) - l + 3)
            Neighbour_Nodes(size(Neighbour_Nodes) - l + 3) = temp
        end do

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8
        r = 0.0_8

        Shape_Functions = 0.0_8

        Shape_Functions_Neighbour = 0.0_8

        do l = 1, Num_Gauss_Points

            do p = 1, size(Nodes)
                k = Properties%Elements(i)%Side_Nodes(j,p)
                Shape_Functions(k,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))
                k_n = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),p)
                Shape_Functions_Neighbour(k_n,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Neighbour_Nodes(p), 1.0_8, xi(l))
                dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
                r(l) = r(l) + Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
            end do

        end do

        F_in = 0.0_8

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

            do b = 1, size(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),:))

                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),b)

                if(Properties%g == 0) F_in(index_1,index_2) = F_in(index_1,index_2) + sum(w*sqrt(dx_dxi**2 + dy_dxi**2)*Shape_Functions(index_1,:)*Shape_Functions_Neighbour(index_2,:))
                if(Properties%g == 1) F_in(index_1,index_2) = F_in(index_1,index_2) + sum(w*sqrt(dx_dxi**2 + dy_dxi**2)*Shape_Functions(index_1,:)*Shape_Functions_Neighbour(index_2,:))

            end do

        end do

    end Function Integrate_Quad_Side_F_in

    subroutine Integrate_Quad_Side_C(Properties,N,i,j,l,F_out,Omega_n)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j, l
        real(kind = 8)      :: Omega_n
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, p
        real(kind = 8), dimension(:,:) :: F_out
        real(kind = 8), dimension(:), allocatable :: xi, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi

        Num_Gauss_Points = (2*N%Degree+2)

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(Nodes(N%Degree+1))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do p = 3, N%Degree+1
            if (p == 3) then
                Nodes(p) = 5
            else
                Nodes(p) = 1 + Nodes(p-1)
            end if
        end do

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8

        Shape_Functions = 0.0_8

        do p = 1, size(Nodes)
            k = Properties%Elements(i)%Side_Nodes(j,p)
            Shape_Functions(k,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))
            dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
            dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
        end do

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

           do b = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(i)%Side_Nodes(j,b)

                F_out(index_1,index_2) = F_out(index_1,index_2) + abs(Omega_n)*w(l)*sqrt(dx_dxi(l)**2 + dy_dxi(l)**2)*Shape_Functions(index_1,l)*Shape_Functions(index_2,l)

            end do

        end do

    end subroutine Integrate_Quad_Side_C

    subroutine Integrate_Quad_Side_F_in_C(Properties,N,i,j,l,F_in,Omega_n)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i, j, l
        real(kind = 8), intent(in) :: Omega_n
        integer             :: Num_Gauss_Points, a, b, k, k_n, p, temp, index_1, index_2

        real(kind = 8), dimension(:,:) :: F_in
        real(kind = 8), dimension(:), allocatable :: xi, w

        integer, dimension(N%Degree+1) :: Nodes, Neighbour_Nodes

        real(kind = 8), dimension(:), allocatable  :: dx_dxi, dy_dxi

        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions, Shape_Functions_Neighbour

        Num_Gauss_Points = (2*N%Degree+2)

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes, Num_Gauss_Points))

        allocate(Shape_Functions_Neighbour(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes, Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do p = 3, N%Degree+1
            if (p == 3) then
                Nodes(p) = 5
            else
                Nodes(p) = 1 + Nodes(p-1)
            end if
        end do

        Neighbour_Nodes = Nodes
        temp = Neighbour_Nodes(1)
        Neighbour_Nodes(1) = Neighbour_Nodes(2)
        Neighbour_Nodes(2) = temp
        do p = 3, 2 + (size(Neighbour_Nodes) - 2)/2
            temp = Neighbour_Nodes(p)
            Neighbour_Nodes(p) = Neighbour_Nodes(size(Neighbour_Nodes) - p + 3)
            Neighbour_Nodes(size(Neighbour_Nodes) - p + 3) = temp
        end do

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8

        Shape_Functions = 0.0_8

        Shape_Functions_Neighbour = 0.0_8

        do p = 1, size(Nodes)
            k = Properties%Elements(i)%Side_Nodes(j,p)
            Shape_Functions(k,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))
            k_n = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),p)
            Shape_Functions_Neighbour(k_n,l) = Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, Neighbour_Nodes(p), 1.0_8, xi(l))
            dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
            dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
        end do

        do a = 1, size(Properties%Elements(i)%Side_Nodes(j,:))

            do b = 1, size(Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),:))
 
                index_1 = Properties%Elements(i)%Side_Nodes(j,a)
                index_2 = Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Side_Nodes(Properties%Elements(i)%Neighbours(j,2),b)

                F_in(index_1,index_2) = F_in(index_1,index_2) + abs(Omega_n)*w(l)*sqrt(dx_dxi(l)**2 + dy_dxi(l)**2)*Shape_Functions(index_1,l)*Shape_Functions_Neighbour(index_2,l)

            end do

        end do

    end subroutine Integrate_Quad_Side_F_in_C
    
end module m_Create_Quadrilateral_Shape_Functions