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

    Function Generate_Shape_Functions(Coordinates, N, a, x) Result(Value)

        type(NTypeD) :: N

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a
        real(kind = 8), dimension(:)              :: Coordinates
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x

        allocate(Shape_Functions(N%Degree+1))

        ! Linear shape functions

        if (N%Degree == 1) then

            Shape_Functions(1) = (x - Coordinates(2)) / (Coordinates(1) - Coordinates(2))
            Shape_Functions(2) = (x - Coordinates(1)) / (Coordinates(2) - Coordinates(1))

        else if (N%Degree == 2) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)))

        else if (N%Degree == 3) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)))

        else if (N%Degree == 4) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(5)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)))
            Shape_Functions(5) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)))

        else if (N%Degree == 5) then

            Shape_Functions(1) = (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)) * (Coordinates(1) - Coordinates(6)))
            Shape_Functions(2) = (x - Coordinates(1)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)) * (Coordinates(2) - Coordinates(6)))
            Shape_Functions(3) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(4)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)) * (Coordinates(3) - Coordinates(6)))
            Shape_Functions(4) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(5)) * (x - Coordinates(6)) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)) * (Coordinates(4) - Coordinates(6)))
            Shape_Functions(5) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(6)) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)) * (Coordinates(5) - Coordinates(6)))
            Shape_Functions(6) = (x - Coordinates(1)) * (x - Coordinates(2)) * (x - Coordinates(3)) * (x - Coordinates(4)) * (x - Coordinates(5)) / ((Coordinates(6) - Coordinates(1)) * (Coordinates(6) - Coordinates(2)) * (Coordinates(6) - Coordinates(3)) * (Coordinates(6) - Coordinates(4)) * (Coordinates(6) - Coordinates(5)))
        
        end if

        Value = Shape_Functions(a)

    end Function Generate_Shape_Functions

    Function Generate_Shape_Functions_Derivative(Coordinates, N, a, x) Result(Value)

        type(NTypeD) :: N

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:)              :: Coordinates
        integer, intent(in)                       :: a
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: x

        allocate(Shape_Functions_Derivative(N%Degree+1))

        ! Linear shape functions

        if (N%Degree == 1) then

            Shape_Functions_Derivative(1) = 1.0 / (Coordinates(1) - Coordinates(2))
            Shape_Functions_Derivative(2) = 1.0 / (Coordinates(2) - Coordinates(1))

        else if (N%Degree == 2) then

            Shape_Functions_Derivative(1) = (2.0 * x - Coordinates(2) - Coordinates(3)) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)))
            Shape_Functions_Derivative(2) = (2.0 * x - Coordinates(1) - Coordinates(3)) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)))
            Shape_Functions_Derivative(3) = (2.0 * x - Coordinates(1) - Coordinates(2)) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)))

        else if (N%Degree == 3) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(4))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2)) + (x - Coordinates(1))*(x - Coordinates(3)) + (x - Coordinates(2))*(x - Coordinates(3))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)))

        else if (N%Degree == 4) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)))
            Shape_Functions_Derivative(5) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)))

        else if (N%Degree == 5) then

            Shape_Functions_Derivative(1) = ((x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(1) - Coordinates(2)) * (Coordinates(1) - Coordinates(3)) * (Coordinates(1) - Coordinates(4)) * (Coordinates(1) - Coordinates(5)) * (Coordinates(1) - Coordinates(6)))
            Shape_Functions_Derivative(2) = ((x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(2) - Coordinates(1)) * (Coordinates(2) - Coordinates(3)) * (Coordinates(2) - Coordinates(4)) * (Coordinates(2) - Coordinates(5)) * (Coordinates(2) - Coordinates(6)))
            Shape_Functions_Derivative(3) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(3) - Coordinates(1)) * (Coordinates(3) - Coordinates(2)) * (Coordinates(3) - Coordinates(4)) * (Coordinates(3) - Coordinates(5)) * (Coordinates(3) - Coordinates(6)))
            Shape_Functions_Derivative(4) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5))*(x - Coordinates(6))) / ((Coordinates(4) - Coordinates(1)) * (Coordinates(4) - Coordinates(2)) * (Coordinates(4) - Coordinates(3)) * (Coordinates(4) - Coordinates(5)) * (Coordinates(4) - Coordinates(6)))
            Shape_Functions_Derivative(5) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(6))) / ((Coordinates(5) - Coordinates(1)) * (Coordinates(5) - Coordinates(2)) * (Coordinates(5) - Coordinates(3)) * (Coordinates(5) - Coordinates(4)) * (Coordinates(5) - Coordinates(6)))
            Shape_Functions_Derivative(6) = ((x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(2))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(1))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5)) + (x - Coordinates(2))*(x - Coordinates(3))*(x - Coordinates(4))*(x - Coordinates(5))) / ((Coordinates(6) - Coordinates(1)) * (Coordinates(6) - Coordinates(2)) * (Coordinates(6) - Coordinates(3)) * (Coordinates(6) - Coordinates(4)) * (Coordinates(6) - Coordinates(5)))
        
        end if

        Value = Shape_Functions_Derivative(a)

    end Function Generate_Shape_Functions_Derivative

    subroutine Calculate_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+1) :: x, w

        real(kind = 8), dimension(N%Degree+1) :: Coordinates

        real(kind = 8), dimension(:,:) :: D_Matrix

        D_Matrix = 0.0_8

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

            Coordinates(1) = -1.0_8

            do j = 2, N%Degree + 1

                Coordinates(j) =  Coordinates(j-1) + 2.0_8/N%Degree

            end do

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+1

                        D_Matrix(index_1,index_2) = D_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions_Derivative(Coordinates, N, index_1, x(j))*Generate_Shape_Functions_Derivative(Coordinates, N, index_2, x(j))

                    end do

                end do

            end do

            D_Matrix = D_Matrix*2.0_8/Properties%Elements(i)%Volume

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(N%Degree+1, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+1

                        D_Matrix(index_1,index_2) = D_Matrix(index_1,index_2) + w(j)*x(j)*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions_Derivative(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do

        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(N%Degree+1, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+1

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
        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+2) :: x, w

        real(kind = 8), dimension(N%Degree+1) :: Coordinates

        real(kind = 8), dimension(:,:) :: Mass_Matrix

        Mass_Matrix = 0.0_8

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(N%Degree+2, x, w)

            Coordinates(1) = -1.0_8

            do j = 2, N%Degree + 1

                Coordinates(j) =  Coordinates(j-1) + 2.0_8/N%Degree

            end do

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+2

                        Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Shape_Functions(Coordinates, N, index_1, x(j))*Generate_Shape_Functions(Coordinates, N, index_2, x(j))

                    end do

                end do

            end do

            Mass_Matrix = Mass_Matrix*Properties%Elements(i)%Volume/2.0_8

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(N%Degree+2, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+2

                        Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*x(j)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_2, x(j))

                    end do

                end do

            end do
        
        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(N%Degree+2, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    do j = 1, N%Degree+2

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
        integer             :: j, index_1

        real(kind = 8), dimension(N%Degree+2) :: x, w

        real(kind = 8), dimension(N%Degree+1) :: Coordinates

        real(kind = 8), dimension(:) :: Source_Vector

        Source_Vector = 0.0_8

        if (Properties%g == 0) then

            call Generate_1D_Quad_Gauss_Points(N%Degree+2, x, w)

            Coordinates(1) = -1.0_8

            do j = 2, N%Degree + 1

                Coordinates(j) =  Coordinates(j-1) + 2.0_8/N%Degree

            end do

            do index_1 = 1, N%Degree+1

                do j = 1, N%Degree+2

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Shape_Functions(Coordinates, N, index_1, x(j))

                end do

            end do

            Source_Vector = Source_Vector*Properties%Elements(i)%Volume/2.0_8

        else if(Properties%g == 1) then

            call Generate_Modified_Gauss_Points(N%Degree+2, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do j = 1, N%Degree+2

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*x(j)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))

                end do

            end do
        
        else if(Properties%g == 2) then

            call Generate_Modified_Gauss_Points(N%Degree+2, x, w, Properties%Elements(i)%Volume, Properties%Elements(i)%Coordinates(1,1)+Properties%Elements(i)%Volume/2.0_8)

            do index_1 = 1, N%Degree+1

                do j = 1, N%Degree+2

                    Source_Vector(index_1) = Source_Vector(index_1) + w(j)*(x(j)**2)*Generate_Shape_Functions(Properties%Elements(i)%Coordinates(:,1), N, index_1, x(j))

                end do

            end do

        end if

    end subroutine Calculate_Source_Vector

end module m_Create_Shape_Functions_D