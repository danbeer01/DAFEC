module m_Create_Shape_Functions_D
!
! Purpose:
! To create shape functions and their derivatives for a line element
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Gauss_Points

implicit none

contains

    subroutine Calculate_Isoparametric_Linear_Coordinates(Degree, Properties)

        type(PropertiesTypeD) :: Properties

        integer, intent(in) :: Degree

        integer :: j

        allocate(Properties%Isoparametric_Coordinates(Degree+1,3))

        Properties%Isoparametric_Coordinates(1,:) = [-1.0_8, 0.0_8, 0.0_8]
        Properties%Isoparametric_Coordinates(2,:) = [1.0_8, 0.0_8, 0.0_8]

        do j = 3, Degree+1
            if (j == 3) then
                Properties%Isoparametric_Coordinates(j,1) = -1.0_8 + 2.0_8/Degree
                Properties%Isoparametric_Coordinates(j,2) = 0.0_8
                Properties%Isoparametric_Coordinates(j,3) = 0.0_8
            else
                Properties%Isoparametric_Coordinates(j,1) = Properties%Isoparametric_Coordinates(j-1,1) + 2.0_8/Degree
                Properties%Isoparametric_Coordinates(j,2) = 0.0_8
                Properties%Isoparametric_Coordinates(j,3) = 0.0_8
            end if
        end do

    end subroutine

    Function Generate_Shape_Functions(Coordinates, N, a, x) Result(Value)

        type(NTypeD) :: N

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a
        real(kind = 8), dimension(:)              :: Coordinates
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x
        integer                                   :: i

        allocate(Shape_Functions(N%Degree+1))

        ! Linear shape functions

        if (N%Degree == 1) then

            Shape_Functions(1) = (x - Coordinates(2)) / (Coordinates(1) - Coordinates(2))
            Shape_Functions(2) = (x - Coordinates(1)) / (Coordinates(2) - Coordinates(1))
            
            Value = Shape_Functions(a)

        else if (N%Degree == 2) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)))
            
            Value = Shape_Functions(a)

        else if (N%Degree == 3) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)))
            
            Value = Shape_Functions(a)

        else if (N%Degree == 4) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(5)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)))
            Shape_Functions(5) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)))
           
            Value = Shape_Functions(a)

        else if (N%Degree == 5) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)) * (Coordinates(1) - Coordinates(6)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)) * (Coordinates(2) - Coordinates(6)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)) * (Coordinates(3) - Coordinates(6)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)) * (Coordinates(4) - Coordinates(6)))
            Shape_Functions(5) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(6)) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)) * (Coordinates(5) - Coordinates(6)))
            Shape_Functions(6) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(6) - Coordinates(1)) * (Coordinates(6) - Coordinates(2)) * (Coordinates(6) - Coordinates(3)) * (Coordinates(6) - Coordinates(4)) * (Coordinates(6) - Coordinates(5)))

            Value = Shape_Functions(a)

        else

            Value = 1.0_8

            do i = 1, N%Degree+1
                if (i /= a) then
                    Value = Value * (x - Coordinates(i)) / (Coordinates(a) - Coordinates(i))
                end if
            end do

        end if

    end Function Generate_Shape_Functions

    Function Generate_Shape_Functions_Derivative(Coordinates, N, a, x) Result(Value)

        type(NTypeD) :: N

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:)              :: Coordinates
        integer, intent(in)                       :: a
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x
        integer                                   :: i, j
        real(kind = 8)                            :: product, product_2, sum

        allocate(Shape_Functions_Derivative(N%Degree+1))

        ! Linear shape functions

        if (N%Degree == 1) then

            Shape_Functions_Derivative(1) = 1.0 / (Coordinates(1) - Coordinates(2))
            Shape_Functions_Derivative(2) = 1.0 / (Coordinates(2) - Coordinates(1))

            Value = Shape_Functions_Derivative(a)

        else if (N%Degree == 2) then

            Shape_Functions_Derivative(1) = (2.0 * x - Coordinates(2) - Coordinates(3)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)))
            Shape_Functions_Derivative(2) = (2.0 * x - Coordinates(1) - Coordinates(3)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)))
            Shape_Functions_Derivative(3) = (2.0 * x - Coordinates(1) - Coordinates(2)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)))

            Value = Shape_Functions_Derivative(a)

        else if (N%Degree == 3) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(4))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(3))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)))

            Value = Shape_Functions_Derivative(a)

        else if (N%Degree == 4) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)))
            Shape_Functions_Derivative(5) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)))

            Value = Shape_Functions_Derivative(a)

        else if (N%Degree == 5) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)) * (Coordinates(1) - Coordinates(6)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)) * (Coordinates(2) - Coordinates(6)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)) * (Coordinates(3) - Coordinates(6)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)) * (Coordinates(4) - Coordinates(6)))
            Shape_Functions_Derivative(5) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6))) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)) * (Coordinates(5) - Coordinates(6)))
            Shape_Functions_Derivative(6) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(6) - Coordinates(1)) * (Coordinates(6) - Coordinates(2)) * (Coordinates(6) - Coordinates(3)) * (Coordinates(6) - Coordinates(4)) * (Coordinates(6) - Coordinates(5)))

            Value = Shape_Functions_Derivative(a)

        else

            product = 1.0_8
            sum = 0.0_8

            do i = 1, N%Degree+1
                if (i /= a) then
                    product = product * (Coordinates(a) - Coordinates(i))
                    product_2 = 1.0_8
                    do j = 1, N%Degree+1
                        if (j /= a .and. j /=i) then
                            product_2 = product_2 * (x - Coordinates(j))
                        end if
                    end do
                    sum = sum + product_2
                end if
            end do

            Value = sum / product

        end if

    end Function Generate_Shape_Functions_Derivative

    subroutine Calculate_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, index_2, Num_gauss_points

        real(kind = 8), dimension(N%Degree+1) :: x, w

        real(kind = 8), dimension(:,:) :: D_Matrix

        D_Matrix = 0.0_8

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points - 1

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(Num_gauss_points, x, w)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        D_Matrix(index_1,index_2) = D_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

            D_Matrix = D_Matrix*2.0_8/Properties%Elements(i)%Volume

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        D_Matrix(index_1,index_2) = D_Matrix(index_1,index_2) + w(j)*x(j)*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        D_Matrix(index_1,index_2) = D_Matrix(index_1,index_2) + w(j)*(x(j)**2)*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

        end if

    end subroutine Calculate_D_Matrix

    subroutine Calculate_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, index_2, Num_gauss_points

        real(kind = 8), dimension(N%Degree+2) :: x, w

        real(kind = 8), dimension(:,:) :: Mass_Matrix

        Mass_Matrix = 0.0_8

        Num_gauss_points = Properties%Elements(i)%Number_of_Gauss_Points

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(Num_gauss_points, x, w)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

            Mass_Matrix = Mass_Matrix*Properties%Elements(i)%Volume/2.0_8

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*x(j)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do
        
        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, Num_gauss_points

                        Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*(x(j)**2)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

        end if

    end subroutine Calculate_Mass_Matrix

    subroutine Calculate_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, Num_gauss_points

        real(kind = 8), dimension(N%Degree+2) :: x, w

        real(kind = 8), dimension(:) :: Source_Vector

        Source_Vector = 0.0_8

        Num_gauss_points = Properties%Elements(i)%Number_of_Gauss_Points

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(Num_gauss_points, x, w)

            do index_1 = 1, N%Degree+1

                do j = 1, Num_gauss_points

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))

                end do

            end do

            Source_Vector = Source_Vector*Properties%Elements(i)%Volume/2.0_8

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do j = 1, Num_gauss_points

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*x(j)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))

                end do

            end do
        
        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(Num_gauss_points, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do j = 1, Num_gauss_points

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*(x(j)**2)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))

                end do

            end do

        end if

    end subroutine Calculate_Source_Vector

end module m_Create_Shape_Functions_D