module m_Create_Shape_Functions
!
! Purpose:
! To create the shape functions for a line element
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

    subroutine Calculate_Isoparametric_Coordinates(Degree, Properties)

        type(PropertiesType) :: Properties

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

    end subroutine Calculate_Isoparametric_Coordinates

    Function Generate_Shape_Functions(Coordinates, N, a, x) Result(Value)

        type(NType) :: N

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a
        real(kind = 8), dimension(:)              :: Coordinates
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x
        integer                                   :: i

        allocate(Shape_Functions(N%Degree+1))

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

    Function Generate_Shape_Functions_Derivative(Coordinates, Degree, a, x) Result(Value)

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:)              :: Coordinates
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x
        integer                                   :: i, j
        real(kind = 8)                            :: product, product_2, sum

        allocate(Shape_Functions_Derivative(Degree+1))

        ! Linear shape functions

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 1.0 / (Coordinates(1) - Coordinates(2))
            Shape_Functions_Derivative(2) = 1.0 / (Coordinates(2) - Coordinates(1))

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            Shape_Functions_Derivative(1) = (2.0 * x - Coordinates(2) - Coordinates(3)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)))
            Shape_Functions_Derivative(2) = (2.0 * x - Coordinates(1) - Coordinates(3)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)))
            Shape_Functions_Derivative(3) = (2.0 * x - Coordinates(1) - Coordinates(2)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)))

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 3) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(4))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(3))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)))

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 4) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)))
            Shape_Functions_Derivative(5) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)))

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 5) then

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

            do i = 1, Degree+1
                if (i /= a) then
                    product = product * (Coordinates(a) - Coordinates(i))
                    product_2 = 1.0_8
                    do j = 1, Degree+1
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

    subroutine Calculate_Streaming_Matrix(Properties,N,i,Streaming_Matrix)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, index_2

        real(kind = 8), dimension(:,:), intent(inout) :: Streaming_Matrix

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8), dimension(N%Degree+1) :: r

        Streaming_Matrix = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1))/2.0_8 + x*(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))/2.0_8
        
        do index_1 = 1, N%Degree+1

            do index_2 = 1, N%Degree+1

                do j = 1, N%Degree+1

                    if (Properties%g == 0) Streaming_Matrix(index_1,index_2) = Streaming_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N%Degree, index_1, x(j))*Properties%Elements(i)%Inverse_Jacobian(j,1,1)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))*Properties%Elements(i)%Det_Jacobian(j)
                    if (Properties%g == 1) Streaming_Matrix(index_1,index_2) = Streaming_Matrix(index_1,index_2) + w(j)*(r(j))*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N%Degree, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))
                    if (Properties%g == 2) Streaming_Matrix(index_1,index_2) = Streaming_Matrix(index_1,index_2) + w(j)*(r(j)**2)*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N%Degree, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))

                end do

            end do

        end do

    end subroutine Calculate_Streaming_Matrix

    subroutine Calculate_Mass_Matrix(Properties,N,i,Mass_Matrix)
        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, Number_Of_Gauss_Points, index_1, index_2

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: x, w, r

        Mass_Matrix = 0.0_8

        if(Properties%g ==0) Number_Of_Gauss_Points = N%Degree+1
        if(Properties%g ==1) Number_Of_Gauss_Points = N%Degree+2
        if(Properties%g ==2) Number_Of_Gauss_Points = N%Degree+2

        allocate(x(Number_Of_Gauss_Points), w(Number_Of_Gauss_Points), r(Number_Of_Gauss_Points))

        call Generate_1D_Quad_Gauss_Points(Number_Of_Gauss_Points, x, w)

        r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1))/2.0_8 + x*abs(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))/2.0_8

        do index_1 = 1, N%Degree+1

            do index_2 = 1, N%Degree+1

                do j = 1, Number_Of_Gauss_Points

                    if (Properties%g == 0) Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))*Properties%Elements(i)%Det_Jacobian(j)
                    if (Properties%g == 1) Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*(r(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))
                    if (Properties%g == 2) Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*(r(j)**2)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x(j))

                end do

            end do

        end do

    end subroutine Calculate_Mass_Matrix

    subroutine Calculate_1D_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8), dimension(N%Degree+1) :: r

        Source_Vector = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1))/2.0_8 + x*(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))/2.0_8

        do index_1 = 1, N%Degree+1

            do j = 1, N%Degree+1

                if (Properties%g == 0) Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))*Properties%Elements(i)%Det_Jacobian(j)
                if (Properties%g == 1) Source_Vector(index_1) = Source_Vector(index_1) + w(j)*(r(j))*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))
                if (Properties%g == 2) Source_Vector(index_1) = Source_Vector(index_1) + w(j)*(r(j)**2)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x(j))

            end do

        end do

    end subroutine Calculate_1D_Source_Vector

end module m_Create_Shape_Functions