module m_Create_Hexahedral_Shape_Functions_D
!
! Purpose:
! To create the shape functions for a hexahedral element
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

    subroutine Calculate_Isoparametric_Hexahedral_Coordinates(Degree, this)

        type(PropertiesTypeD), intent(inout) :: this

        real(kind = 8)                              :: x, y
        integer, intent(in)                         :: Degree
        integer :: i, j, k, Coord_index

        allocate(this%Isoparametric_Coordinates((Degree+1)*(Degree+1)*(Degree+1),3))

        Coord_index = 1

        k = 0

        if (Degree == 1) then

            ! Corner nodes
            this%Isoparametric_Coordinates(1,:) = [1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(2,:) = [-1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(3,:) = [-1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(4,:) = [1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(5,:) = [1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(6,:) = [-1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(7,:) = [-1.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(8,:) = [1.0_8, -1.0_8, -1.0_8]

        else if (Degree == 2) then
            
            x = 1.0_8 - 2.0_8*k/real(Degree,kind=8)
            y = 1.0_8 - 2.0_8*k/real(Degree,kind=8)
            i = 1
            j = 1

            ! Corner nodes
            this%Isoparametric_Coordinates(1,:) = [1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(2,:) = [-1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(3,:) = [-1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(4,:) = [1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(5,:) = [1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(6,:) = [-1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(7,:) = [-1.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(8,:) = [1.0_8, -1.0_8, -1.0_8]

            ! Upper face nodes
            this%Isoparametric_Coordinates(9,:) = [0.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(10,:) = [-1.0_8, 0.0_8, 1.0_8]
            this%Isoparametric_Coordinates(11,:) = [0.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(12,:) = [1.0_8, 0.0_8, 1.0_8]

            ! Lower face nodes
            this%Isoparametric_Coordinates(13,:) = [0.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(14,:) = [-1.0_8, 0.0_8, -1.0_8]
            this%Isoparametric_Coordinates(15,:) = [0.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(16,:) = [1.0_8, 0.0_8, -1.0_8]

            ! Middle corner nodes
            this%Isoparametric_Coordinates(17,:) = [1.0_8, 1.0_8, 0.0_8]
            this%Isoparametric_Coordinates(18,:) = [-1.0_8, 1.0_8, 0.0_8]
            this%Isoparametric_Coordinates(19,:) = [-1.0_8, -1.0_8, 0.0_8]
            this%Isoparametric_Coordinates(20,:) = [1.0_8, -1.0_8, 0.0_8]

            ! Middle face nodes
            this%Isoparametric_Coordinates(21,:) = [1.0_8, 0.0_8, 0.0_8]
            this%Isoparametric_Coordinates(22,:) = [-1.0_8, 0.0_8, 0.0_8]
            this%Isoparametric_Coordinates(23,:) = [0.0_8, 1.0_8, 0.0_8]
            this%Isoparametric_Coordinates(24,:) = [0.0_8, -1.0_8, 0.0_8]
            this%Isoparametric_Coordinates(25,:) = [0.0_8, 0.0_8, 1.0_8]
            this%Isoparametric_Coordinates(26,:) = [0.0_8, 0.0_8, -1.0_8]

            this%Isoparametric_Coordinates(27,:) = [0.0_8, 0.0_8, 0.0_8]

        else if (Degree == 3) then

            ! Corner nodes
            this%Isoparametric_Coordinates(1,:) = [1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(2,:) = [-1.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(3,:) = [-1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(4,:) = [1.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(5,:) = [1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(6,:) = [-1.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(7,:) = [-1.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(8,:) = [1.0_8, -1.0_8, -1.0_8]

            ! Upper face side nodes
            this%Isoparametric_Coordinates(9,:) = [1.0_8/3.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(10,:) = [-1.0_8/3.0_8, 1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(11,:) = [-1.0_8, 1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(12,:) = [-1.0_8, -1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(13,:) = [-1.0_8/3.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(14,:) = [1.0_8/3.0_8, -1.0_8, 1.0_8]
            this%Isoparametric_Coordinates(15,:) = [1.0_8, -1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(16,:) = [1.0_8, 1.0_8/3.0_8, 1.0_8]

            ! Lower face side nodes
            this%Isoparametric_Coordinates(17,:) = [1.0_8/3.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(18,:) = [-1.0_8/3.0_8, 1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(19,:) = [-1.0_8, 1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(20,:) = [-1.0_8, -1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(21,:) = [-1.0_8/3.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(22,:) = [1.0_8/3.0_8, -1.0_8, -1.0_8]
            this%Isoparametric_Coordinates(23,:) = [1.0_8, -1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(24,:) = [1.0_8, 1.0_8/3.0_8, -1.0_8]

            ! Middle side nodes
            this%Isoparametric_Coordinates(25,:) = [1.0_8, 1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(26,:) = [1.0_8, 1.0_8, -1.0_8/3.0_8]

            this%Isoparametric_Coordinates(27,:) = [-1.0_8, 1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(28,:) = [-1.0_8, 1.0_8, -1.0_8/3.0_8]

            this%Isoparametric_Coordinates(29,:) = [-1.0_8, -1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(30,:) = [-1.0_8, -1.0_8, -1.0_8/3.0_8]

            this%Isoparametric_Coordinates(31,:) = [1.0_8, -1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(32,:) = [1.0_8, -1.0_8, -1.0_8/3.0_8]

            ! Upper face nodes
            this%Isoparametric_Coordinates(33,:) = [1.0_8/3.0_8, 1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(34,:) = [-1.0_8/3.0_8, 1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(35,:) = [-1.0_8/3.0_8, -1.0_8/3.0_8, 1.0_8]
            this%Isoparametric_Coordinates(36,:) = [1.0_8/3.0_8, -1.0_8/3.0_8, 1.0_8]

            ! Lower face nodes
            this%Isoparametric_Coordinates(37,:) = [1.0_8/3.0_8, 1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(38,:) = [-1.0_8/3.0_8, 1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(39,:) = [-1.0_8/3.0_8, -1.0_8/3.0_8, -1.0_8]
            this%Isoparametric_Coordinates(40,:) = [1.0_8/3.0_8, -1.0_8/3.0_8, -1.0_8]

            ! Middle faces nodes
            this%Isoparametric_Coordinates(41,:) = [1.0_8/3.0_8, 1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(42,:) = [-1.0_8/3.0_8, 1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(43,:) = [-1.0_8/3.0_8, 1.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(44,:) = [1.0_8/3.0_8, 1.0_8, -1.0_8/3.0_8]

            this%Isoparametric_Coordinates(45,:) = [-1.0_8, 1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(46,:) = [-1.0_8, -1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(47,:) = [-1.0_8, -1.0_8/3.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(48,:) = [-1.0_8, 1.0_8/3.0_8, -1.0_8/3.0_8]

            ! 3,4,7,8
            this%Isoparametric_Coordinates(49,:) = [-1.0_8/3.0_8, -1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(50,:) = [1.0_8/3.0_8, -1.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(51,:) = [1.0_8/3.0_8, -1.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(52,:) = [-1.0_8/3.0_8, -1.0_8, -1.0_8/3.0_8]

            ! 1,4,5,8
            this%Isoparametric_Coordinates(53,:) = [1.0_8, -1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(54,:) = [1.0_8, 1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(55,:) = [1.0_8, 1.0_8/3.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(56,:) = [1.0_8, -1.0_8/3.0_8, -1.0_8/3.0_8]

            ! Middle nodes
            this%Isoparametric_Coordinates(57,:) = [1.0_8/3.0_8, 1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(58,:) = [-1.0_8/3.0_8, 1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(59,:) = [-1.0_8/3.0_8, -1.0_8/3.0_8, 1.0_8/3.0_8]
            this%Isoparametric_Coordinates(60,:) = [1.0_8/3.0_8, -1.0_8/3.0_8, 1.0_8/3.0_8]

            this%Isoparametric_Coordinates(61,:) = [1.0_8/3.0_8, 1.0_8/3.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(62,:) = [-1.0_8/3.0_8, 1.0_8/3.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(63,:) = [-1.0_8/3.0_8, -1.0_8/3.0_8, -1.0_8/3.0_8]
            this%Isoparametric_Coordinates(64,:) = [1.0_8/3.0_8, -1.0_8/3.0_8, -1.0_8/3.0_8]

        end if

    end subroutine Calculate_Isoparametric_Hexahedral_Coordinates

    Function Generate_Hexahedral_Shape_Functions(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: eta, xi, zeta

        integer :: i

        Value = 0.0_8

        allocate(Shape_Functions((Degree+1)*(Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions(1) = 0.125_8*(1.0_8 + xi)*(1.0_8 + eta)*(1.0_8 + zeta)
            Shape_Functions(2) = 0.125_8*(1.0_8 - xi)*(1.0_8 + eta)*(1.0_8 + zeta)
            Shape_Functions(3) = 0.125_8*(1.0_8 - xi)*(1.0_8 - eta)*(1.0_8 + zeta)
            Shape_Functions(4) = 0.125_8*(1.0_8 + xi)*(1.0_8 - eta)*(1.0_8 + zeta)
            Shape_Functions(5) = 0.125_8*(1.0_8 + xi)*(1.0_8 + eta)*(1.0_8 - zeta)
            Shape_Functions(6) = 0.125_8*(1.0_8 - xi)*(1.0_8 + eta)*(1.0_8 - zeta)
            Shape_Functions(7) = 0.125_8*(1.0_8 - xi)*(1.0_8 - eta)*(1.0_8 - zeta)
            Shape_Functions(8) = 0.125_8*(1.0_8 + xi)*(1.0_8 - eta)*(1.0_8 - zeta)

            Value = Shape_Functions(a)

        else

            Value = 1.0_8

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(xi - this%Isoparametric_Coordinates(i,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(eta - this%Isoparametric_Coordinates(i,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        Value = Value*(zeta - this%Isoparametric_Coordinates(i,3))/(this%Isoparametric_Coordinates(a,3) - this%Isoparametric_Coordinates(i,3))

                    end if

                end if

            end do

        end if

    end Function Generate_Hexahedral_Shape_Functions

    Function Generate_Shape_Functions_Derivative_xi(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, temp
        real(kind = 8), intent(in)                :: xi, eta, zeta

        integer :: i, j

        Value = 0.0_8

        allocate(Shape_Functions_Derivative((Degree+1)*(Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 0.125_8*(1.0_8 + eta)*(1.0_8 + zeta)
            Shape_Functions_Derivative(2) = -0.125_8*(1.0_8 + eta)*(1.0_8 + zeta)
            Shape_Functions_Derivative(3) = -0.125_8*(1.0_8 - eta)*(1.0_8 + zeta)
            Shape_Functions_Derivative(4) = 0.125_8*(1.0_8 - eta)*(1.0_8 + zeta)
            Shape_Functions_Derivative(5) = 0.125_8*(1.0_8 + eta)*(1.0_8 - zeta)
            Shape_Functions_Derivative(6) = -0.125_8*(1.0_8 + eta)*(1.0_8 - zeta)
            Shape_Functions_Derivative(7) = -0.125_8*(1.0_8 - eta)*(1.0_8 - zeta)
            Shape_Functions_Derivative(8) = 0.125_8*(1.0_8 - eta)*(1.0_8 - zeta)

            Value = Shape_Functions_Derivative(a)

        else

            Value = 0.0_8

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        temp = 1.0_8

                        do j = 1,(Degree+1)*(Degree+1)*(Degree+1)

                            if (j /= a) then

                                if (ABS(this%Isoparametric_Coordinates(j,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(j,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

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

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(eta - this%Isoparametric_Coordinates(i,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        Value = Value*(zeta - this%Isoparametric_Coordinates(i,3))/(this%Isoparametric_Coordinates(a,3) - this%Isoparametric_Coordinates(i,3))

                    end if

                end if

            end do

        end if

    end Function Generate_Shape_Functions_Derivative_xi

    Function Generate_Shape_Functions_Derivative_eta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, temp
        real(kind = 8), intent(in)                :: xi, eta, zeta

        integer :: i, j

        Value = 0.0_8

        allocate(Shape_Functions_Derivative((Degree+1)*(Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 0.125_8*(1.0_8 + xi)*(1.0_8 + zeta)
            Shape_Functions_Derivative(2) = 0.125_8*(1.0_8 - xi)*(1.0_8 + zeta)
            Shape_Functions_Derivative(3) = -0.125_8*(1.0_8 - xi)*(1.0_8 + zeta)
            Shape_Functions_Derivative(4) = -0.125_8*(1.0_8 + xi)*(1.0_8 + zeta)
            Shape_Functions_Derivative(5) = 0.125_8*(1.0_8 + xi)*(1.0_8 - zeta)
            Shape_Functions_Derivative(6) = 0.125_8*(1.0_8 - xi)*(1.0_8 - zeta)
            Shape_Functions_Derivative(7) = -0.125_8*(1.0_8 - xi)*(1.0_8 - zeta)
            Shape_Functions_Derivative(8) = -0.125_8*(1.0_8 + xi)*(1.0_8 - zeta)

            Value = Shape_Functions_Derivative(a)

        else

            Value = 0.0_8

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        temp = 1.0_8

                        do j = 1,(Degree+1)*(Degree+1)*(Degree+1)

                            if (j /= a) then

                                if (ABS(this%Isoparametric_Coordinates(j,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(j,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

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

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(xi - this%Isoparametric_Coordinates(i,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        Value = Value*(zeta - this%Isoparametric_Coordinates(i,3))/(this%Isoparametric_Coordinates(a,3) - this%Isoparametric_Coordinates(i,3))

                    end if

                end if

            end do

        end if

    end Function Generate_Shape_Functions_Derivative_eta

    Function Generate_Shape_Functions_Derivative_zeta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value, temp
        real(kind = 8), intent(in)                :: xi, eta, zeta

        integer :: i, j

        Value = 0.0_8

        allocate(Shape_Functions_Derivative((Degree+1)*(Degree+1)*(Degree+1)))

        if (Degree == 1) then

            Shape_Functions_Derivative(1) = 0.125_8*(1.0_8 + xi)*(1.0_8 + eta)
            Shape_Functions_Derivative(2) = 0.125_8*(1.0_8 - xi)*(1.0_8 + eta)
            Shape_Functions_Derivative(3) = 0.125_8*(1.0_8 - xi)*(1.0_8 - eta)
            Shape_Functions_Derivative(4) = 0.125_8*(1.0_8 + xi)*(1.0_8 - eta)
            Shape_Functions_Derivative(5) = -0.125_8*(1.0_8 + xi)*(1.0_8 + eta)
            Shape_Functions_Derivative(6) = -0.125_8*(1.0_8 - xi)*(1.0_8 + eta)
            Shape_Functions_Derivative(7) = -0.125_8*(1.0_8 - xi)*(1.0_8 - eta)
            Shape_Functions_Derivative(8) = -0.125_8*(1.0_8 + xi)*(1.0_8 - eta)

            Value = Shape_Functions_Derivative(a)

        else

            Value = 0.0_8

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                        temp = 1.0_8

                        do j = 1,(Degree+1)*(Degree+1)*(Degree+1)

                            if (j /= a) then

                                if (ABS(this%Isoparametric_Coordinates(j,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(j,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6) then

                                    if (j /= i) then

                                        temp = temp*(zeta - this%Isoparametric_Coordinates(j,3))/(this%Isoparametric_Coordinates(a,3) - this%Isoparametric_Coordinates(j,3))

                                    end if

                                end if

                            end if

                        end do

                        Value = Value + temp/(this%Isoparametric_Coordinates(a,3) - this%Isoparametric_Coordinates(i,3))

                    end if

                end if

            end do

            do i = 1,(Degree+1)*(Degree+1)*(Degree+1)

                if (i /= a) then

                    if (ABS(this%Isoparametric_Coordinates(i,2) - this%Isoparametric_Coordinates(a,2)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(xi - this%Isoparametric_Coordinates(i,1))/(this%Isoparametric_Coordinates(a,1) - this%Isoparametric_Coordinates(i,1))

                    else if (ABS(this%Isoparametric_Coordinates(i,1) - this%Isoparametric_Coordinates(a,1)) < 1e-6 .and. ABS(this%Isoparametric_Coordinates(i,3) - this%Isoparametric_Coordinates(a,3)) < 1e-6) then

                        Value = Value*(eta - this%Isoparametric_Coordinates(i,2))/(this%Isoparametric_Coordinates(a,2) - this%Isoparametric_Coordinates(i,2))

                    end if

                end if

            end do

        end if

    end Function Generate_Shape_Functions_Derivative_zeta

    subroutine Calculate_Hex_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:) :: D_Matrix

        Num_Gauss_Points = (N%Degree+1)**3

        allocate(dSFMatT(Num_Gauss_Points,Properties%Elements(i)%Number_of_Nodes,3), dSFMat(Num_Gauss_Points,3,Properties%Elements(i)%Number_of_Nodes))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_3D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        D_Matrix = 0.0_8

        do k = 1, Properties%Elements(i)%Number_of_Nodes

            do j = 1, Num_Gauss_Points

                dSFMat(j,1,k) = Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,2,k) = Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,3,k) = Generate_Shape_Functions_Derivative_zeta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))

                dSFMatT(j,k,1) = dSFMat(j,1,k)
                dSFMatT(j,k,2) = dSFMat(j,2,k)
                dSFMatT(j,k,3) = dSFMat(j,3,k)

            end do

        end do

        do j = 1, Num_Gauss_Points

            D_Matrix = D_Matrix + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(matmul(dSFMatT(j,:,:), Properties%Elements(i)%Inverse_Jacobian(j,:,:)), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do
        
    end subroutine Calculate_Hex_D_Matrix

    subroutine Calculate_Hex_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, index_2, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_3D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Mass_Matrix = 0.0_8

        do index_1 = 1, Num_Nodes

            do index_2 = 1, Num_Nodes

                do j = 1, Num_Gauss_Points

                    Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Hexahedral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*Generate_Hexahedral_Shape_Functions(Properties, N%Degree, index_2, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

                end do

            end do

        end do

    end subroutine Calculate_Hex_Mass_Matrix

    subroutine Calculate_Hex_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Nodes

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_3D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Source_Vector = 0.0_8

        do index_1 = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Hexahedral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

            end do

        end do

    end subroutine Calculate_Hex_Source_Vector

    Function Hex_Face_Boundary(Properties,N,i,side_nodes) Result(B_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i

        integer, dimension(:) :: side_nodes

        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p

        real(kind = 8), dimension(:,:), allocatable :: B_Matrix
        real(kind = 8), dimension(:), allocatable :: xi, eta, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi, dz_dxi, dx_deta, dy_deta, dz_deta

        allocate(B_Matrix(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        Num_Gauss_Points = (N%Degree+1)**2

        allocate(Shape_Functions(Properties%Elements(i)%Number_of_Nodes,Num_Gauss_Points))

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points), dx_deta(Num_Gauss_Points), dy_deta(Num_Gauss_Points), dz_dxi(Num_Gauss_Points), dz_deta(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(Nodes((N%Degree+1)**2))

        call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        Nodes(3) = 3
        Nodes(4) = 4
        if (N%Degree == 2) then
            Nodes(5:9) = [9,10,11,12,25]
        else if (N%Degree == 3) then
            Nodes(5:16) = [9,10,11,12,13,14,15,16,33,34,35,36]
        end if

        if (N%Degree == 1) then
            if (side_nodes(1) == 2 .and. side_nodes(2) == 3 .and. side_nodes(3) == 6 .and. side_nodes(4) == 7) side_nodes = (/2,6,7,3/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 4 .and. side_nodes(3) == 5 .and. side_nodes(4) == 8) side_nodes = (/1,5,8,4/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 2 .and. side_nodes(3) == 5 .and. side_nodes(4) == 6) side_nodes = [1,2,6,5]
            if (side_nodes(1) == 3 .and. side_nodes(2) == 4 .and. side_nodes(3) == 7 .and. side_nodes(4) == 8) side_nodes = [4,3,7,8]
        else if (N%Degree == 2) then
            if (side_nodes(1) == 2 .and. side_nodes(2) == 3 .and. side_nodes(3) == 6 .and. side_nodes(4) == 7) side_nodes = (/2,6,7,3,18,14,19,10,22/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 4 .and. side_nodes(3) == 5 .and. side_nodes(4) == 8) side_nodes = (/1,5,8,4,17,16,20,12,21/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 2 .and. side_nodes(3) == 5 .and. side_nodes(4) == 6) side_nodes = (/1,2,6,5,9,18,13,17,23/)
            if (side_nodes(1) == 3 .and. side_nodes(2) == 4 .and. side_nodes(3) == 7 .and. side_nodes(4) == 8) side_nodes = (/4,3,7,8,11,19,15,20,24/)
        else if (N%Degree == 3) then
            if (side_nodes(1) == 2 .and. side_nodes(2) == 3 .and. side_nodes(3) == 6 .and. side_nodes(4) == 7) side_nodes = (/2,6,7,3,27,28,19,20,30,29,12,11,45,48,47,46/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 4 .and. side_nodes(3) == 5 .and. side_nodes(4) == 8) side_nodes = (/1,5,8,4,25,26,24,23,32,31,15,16,54,55,56,53/)
            if (side_nodes(1) == 1 .and. side_nodes(2) == 2 .and. side_nodes(3) == 5 .and. side_nodes(4) == 6) side_nodes = (/1,2,6,5,9,10,27,28,18,17,26,25,41,42,43,44/)
            if (side_nodes(1) == 3 .and. side_nodes(2) == 4 .and. side_nodes(3) == 7 .and. side_nodes(4) == 8) side_nodes = (/4,3,7,8,14,13,29,30,21,22,32,31,50,49,52,51/)
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
                Shape_Functions(k,l) = Generate_Hexahedral_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)
                dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                dz_dxi(l) = dz_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,3)
                dx_deta(l) = dx_deta(l) + Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,1)
                dy_deta(l) = dy_deta(l) + Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,2)
                dz_deta(l) = dz_deta(l) + Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 1.0_8)*Properties%Elements(i)%Coordinates(k,3)
            end do

        end do

        do a = 1, size(side_nodes)

           do b = 1, size(side_nodes)

                index_1 = side_nodes(a)
                index_2 = side_nodes(b)

                B_Matrix(index_1,index_2) = B_Matrix(index_1,index_2) + 0.5_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                
            end do

        end do

    end Function Hex_Face_Boundary

end module