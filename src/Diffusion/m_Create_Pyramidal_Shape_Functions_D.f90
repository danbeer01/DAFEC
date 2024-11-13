module m_Create_Pyramidal_Shape_Functions_D
!
! Purpose:
! To create the shape functions for a Pyramidal element
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Create_Tetrahedral_Shape_Functions_D
use m_Gauss_Points

implicit none

contains

    Function Generate_Pyramidal_Shape_Functions(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(inout) :: this

        real(kind = 8), dimension(:), allocatable   :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: B_Functions
        real(kind = 8), dimension(:,:), allocatable :: VDM_Inverse, VDM
        integer, intent(in)                         :: a, Degree
        real(kind = 8)                              :: Value
        real(kind = 8)                              :: eta, xi, zeta

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions(5))

            Shape_Functions(1) = 0.25_8*(1.0_8 - xi + eta - zeta - (xi*eta)/(1.0_8 - zeta))
            Shape_Functions(4) = 0.25_8*(1.0_8 + xi + eta - zeta + (xi*eta)/(1.0_8 - zeta))
            Shape_Functions(3) = 0.25_8*(1.0_8 + xi - eta - zeta - (xi*eta)/(1.0_8 - zeta))
            Shape_Functions(2) = 0.25_8*(1.0_8 - xi - eta - zeta + (xi*eta)/(1.0_8 - zeta))
            Shape_Functions(5) = zeta

            Value = Shape_Functions(a)

        else if (Degree == 2) then

            allocate(Shape_Functions(14))
            allocate(B_Functions(14))
            allocate(VDM_Inverse(14,14))
            allocate(VDM(14,14))

            VDM_Inverse(1,:) = (/ -0.020833333334583412_8, -0.027777777777777762_8, 0.027777777777777797_8, -0.030092592590509261_8, -0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555566_8, -0.055555555555555559_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(2,:) = (/ -0.020833333334583315_8, 0.027777777777777762_8, 0.027777777777777762_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, -0.055555555555555566_8, -0.055555555555555566_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(3,:) = (/ -0.020833333334583329_8, 0.027777777777777797_8, -0.027777777777777762_8, -0.030092592590509268_8, -0.16666666666666669_8, 0.018518518521851862_8, -0.055555555555555559_8, 0.055555555555555566_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(4,:) = (/ -0.020833333334583315_8, -0.027777777777777797_8, -0.027777777777777797_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555559_8, 0.055555555555555559_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(5,:) = (/ -0.050000000015000032_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.08333333358333331_8, 0.0000000000000000_8, 0.13333333337333336_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8 /)
            VDM_Inverse(6,:) = (/ 0.06666666666666652_8, 0.0_8, 0.27777777777777779_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.0_8, -0.055555555555555532_8, 0.0_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, -0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(7,:) = (/ 0.06666666666666652_8, 0.27777777777777779_8, 0.0000000000000000_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, -0.055555555555555539_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, -0.0000000000000000_8, -0.33333333333333331_8 /)
            VDM_Inverse(8,:) = (/ 0.06666666666666652_8, 0.0000000000000000_8, -0.27777777777777779_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, 0.0000000000000000_8, 0.055555555555555539_8, -0.0000000000000000_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, 0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(9,:) = (/ 0.06666666666666652_8, -0.27777777777777779_8, 0.0_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.055555555555555532_8, 0.0_8, 0.0_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, 0.0000000000000000_8, 0.33333333333333331_8 /)
            VDM_Inverse(10,:) = (/ 0.15000000000500011_8, -0.16666666666666663_8, 0.16666666666666663_8, 0.08333333332500036_8, -0.16666666666666663_8, -0.06666666668000028_8, -0.16666666666666663_8, 0.16666666666666666_8, -0.16666666666666663_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(11,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(12,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, -0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(13,:) = (/ 0.15000000000499999_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(14,:) = (/ 0.26666666666666661_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.14814814814814814_8, 0.0000000000000000_8, 0.029629629629629631_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.44444444444444442_8, -0.44444444444444442_8, 0.44444444444444442_8, 0.0000000000000000_8, 0.0000000000000000_8 /)

            B_Functions(1) = 1.0_8
            B_Functions(2) = xi
            B_Functions(3) = eta
            B_Functions(4) = 4.0_8*zeta - 1.0_8
            B_Functions(5) = xi*eta/(1.0_8 - zeta)
            B_Functions(6) = 15.0_8*zeta**2 - 10.0_8*zeta + 1.0_8
            B_Functions(7) = xi*(6.0_8*zeta - 1.0_8)
            B_Functions(8) = eta*(6.0_8*zeta - 1.0_8)
            B_Functions(9) = xi*eta*(6.0_8*zeta - 1.0_8)/(1.0_8 - zeta)
            B_Functions(10) = 0.5_8*(3.0_8*xi**2 - (1.0_8 - zeta)**2)
            B_Functions(11) = 0.5_8*(3.0_8*eta**2 - (1.0_8 - zeta)**2)
            B_Functions(12) = 0.25_8*(3.0_8*xi**2/(1.0_8 - zeta) - (1.0_8 - zeta))*(3.0_8*eta**2/(1.0_8 - zeta) - (1.0_8 - zeta))
            B_Functions(13) = 0.5_8*(3.0_8*xi**2 - (1.0_8 - zeta)**2)*eta/(1.0_8 - zeta)
            B_Functions(14) = 0.5_8*(3.0_8*eta**2 - (1.0_8 - zeta)**2)*xi/(1.0_8 - zeta)

            Shape_Functions(:) = matmul(VDM_Inverse, B_Functions)

            Value = Shape_Functions(a)

        end if

    end Function Generate_Pyramidal_Shape_Functions

    Function Generate_Pyramidal_Shape_Functions_Derivative_xi(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:), allocatable :: B_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta
        real(kind = 8), dimension(:,:), allocatable :: VDM_Inverse

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(5))

            Shape_Functions_Derivative(1) = 0.25_8*(-1.0_8 - eta/(1.0_8 - zeta))
            Shape_Functions_Derivative(4) = 0.25_8*(1.0_8 + eta/(1.0_8 - zeta))
            Shape_Functions_Derivative(3) = 0.25_8*(1.0_8 - eta/(1.0_8 - zeta))
            Shape_Functions_Derivative(2) = 0.25_8*(-1.0_8 + eta/(1.0_8 - zeta))
            Shape_Functions_Derivative(5) = 0.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(14))
            allocate(B_Functions_Derivative(14))
            allocate(VDM_Inverse(14,14))

            VDM_Inverse(1,:) = (/ -0.020833333334583412_8, -0.027777777777777762_8, 0.027777777777777797_8, -0.030092592590509261_8, -0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555566_8, -0.055555555555555559_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(2,:) = (/ -0.020833333334583315_8, 0.027777777777777762_8, 0.027777777777777762_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, -0.055555555555555566_8, -0.055555555555555566_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(3,:) = (/ -0.020833333334583329_8, 0.027777777777777797_8, -0.027777777777777762_8, -0.030092592590509268_8, -0.16666666666666669_8, 0.018518518521851862_8, -0.055555555555555559_8, 0.055555555555555566_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(4,:) = (/ -0.020833333334583315_8, -0.027777777777777797_8, -0.027777777777777797_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555559_8, 0.055555555555555559_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(5,:) = (/ -0.050000000015000032_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.08333333358333331_8, 0.0000000000000000_8, 0.13333333337333336_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8 /)
            VDM_Inverse(6,:) = (/ 0.06666666666666652_8, 0.0_8, 0.27777777777777779_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.0_8, -0.055555555555555532_8, 0.0_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, -0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(7,:) = (/ 0.06666666666666652_8, 0.27777777777777779_8, 0.0000000000000000_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, -0.055555555555555539_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, -0.0000000000000000_8, -0.33333333333333331_8 /)
            VDM_Inverse(8,:) = (/ 0.06666666666666652_8, 0.0000000000000000_8, -0.27777777777777779_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, 0.0000000000000000_8, 0.055555555555555539_8, -0.0000000000000000_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, 0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(9,:) = (/ 0.06666666666666652_8, -0.27777777777777779_8, 0.0_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.055555555555555532_8, 0.0_8, 0.0_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, 0.0000000000000000_8, 0.33333333333333331_8 /)
            VDM_Inverse(10,:) = (/ 0.15000000000500011_8, -0.16666666666666663_8, 0.16666666666666663_8, 0.08333333332500036_8, -0.16666666666666663_8, -0.06666666668000028_8, -0.16666666666666663_8, 0.16666666666666666_8, -0.16666666666666663_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(11,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(12,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, -0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(13,:) = (/ 0.15000000000499999_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(14,:) = (/ 0.26666666666666661_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.14814814814814814_8, 0.0000000000000000_8, 0.029629629629629631_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.44444444444444442_8, -0.44444444444444442_8, 0.44444444444444442_8, 0.0000000000000000_8, 0.0000000000000000_8 /)

            B_Functions_Derivative(1) = 0.0_8
            B_Functions_Derivative(2) = 1.0_8
            B_Functions_Derivative(3) = 0.0_8
            B_Functions_Derivative(4) = 0.0_8
            B_Functions_Derivative(5) = eta/(1.0_8 - zeta)
            B_Functions_Derivative(6) = 0.0_8
            B_Functions_Derivative(7) = (6.0_8*zeta - 1.0_8)
            B_Functions_Derivative(8) = 0.0_8
            B_Functions_Derivative(9) = eta*(6.0_8*zeta - 1.0_8)/(1.0_8 - zeta)
            B_Functions_Derivative(10) = 3.0_8*xi
            B_Functions_Derivative(11) = 0.0_8
            B_Functions_Derivative(12) = 0.5_8*(3.0_8*xi/(1.0_8 - zeta))*(3.0_8*eta**2/(1.0_8 - zeta) - (1.0_8 - zeta))
            B_Functions_Derivative(13) = (3.0_8*xi)*eta/(1.0_8 - zeta)
            B_Functions_Derivative(14) = 0.5_8*(3.0_8*eta**2 - (1.0_8 - zeta)**2)*1.0_8/(1.0_8 - zeta)

            Shape_Functions_Derivative(:) = matmul(VDM_Inverse, B_Functions_Derivative)

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Pyramidal_Shape_Functions_Derivative_xi

    Function Generate_Pyramidal_Shape_Functions_Derivative_eta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:), allocatable :: B_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta
        real(kind = 8), dimension(:,:), allocatable :: VDM_Inverse

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(5))

            Shape_Functions_Derivative(1) = 0.25_8*(1.0_8 - xi/(1.0_8 - zeta))
            Shape_Functions_Derivative(4) = 0.25_8*(1.0_8 + xi/(1.0_8 - zeta))
            Shape_Functions_Derivative(3) = 0.25_8*(-1.0_8 - xi/(1.0_8 - zeta))
            Shape_Functions_Derivative(2) = 0.25_8*(-1.0_8 + xi/(1.0_8 - zeta))
            Shape_Functions_Derivative(5) = 0.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(14))
            allocate(B_Functions_Derivative(14))
            allocate(VDM_Inverse(14,14))

            VDM_Inverse(1,:) = (/ -0.020833333334583412_8, -0.027777777777777762_8, 0.027777777777777797_8, -0.030092592590509261_8, -0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555566_8, -0.055555555555555559_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(2,:) = (/ -0.020833333334583315_8, 0.027777777777777762_8, 0.027777777777777762_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, -0.055555555555555566_8, -0.055555555555555566_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(3,:) = (/ -0.020833333334583329_8, 0.027777777777777797_8, -0.027777777777777762_8, -0.030092592590509268_8, -0.16666666666666669_8, 0.018518518521851862_8, -0.055555555555555559_8, 0.055555555555555566_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(4,:) = (/ -0.020833333334583315_8, -0.027777777777777797_8, -0.027777777777777797_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555559_8, 0.055555555555555559_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(5,:) = (/ -0.050000000015000032_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.08333333358333331_8, 0.0000000000000000_8, 0.13333333337333336_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8 /)
            VDM_Inverse(6,:) = (/ 0.06666666666666652_8, 0.0_8, 0.27777777777777779_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.0_8, -0.055555555555555532_8, 0.0_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, -0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(7,:) = (/ 0.06666666666666652_8, 0.27777777777777779_8, 0.0000000000000000_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, -0.055555555555555539_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, -0.0000000000000000_8, -0.33333333333333331_8 /)
            VDM_Inverse(8,:) = (/ 0.06666666666666652_8, 0.0000000000000000_8, -0.27777777777777779_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, 0.0000000000000000_8, 0.055555555555555539_8, -0.0000000000000000_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, 0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(9,:) = (/ 0.06666666666666652_8, -0.27777777777777779_8, 0.0_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.055555555555555532_8, 0.0_8, 0.0_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, 0.0000000000000000_8, 0.33333333333333331_8 /)
            VDM_Inverse(10,:) = (/ 0.15000000000500011_8, -0.16666666666666663_8, 0.16666666666666663_8, 0.08333333332500036_8, -0.16666666666666663_8, -0.06666666668000028_8, -0.16666666666666663_8, 0.16666666666666666_8, -0.16666666666666663_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(11,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(12,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, -0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(13,:) = (/ 0.15000000000499999_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(14,:) = (/ 0.26666666666666661_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.14814814814814814_8, 0.0000000000000000_8, 0.029629629629629631_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.44444444444444442_8, -0.44444444444444442_8, 0.44444444444444442_8, 0.0000000000000000_8, 0.0000000000000000_8 /)

            B_Functions_Derivative(1) = 0.0_8
            B_Functions_Derivative(2) = 0.0_8
            B_Functions_Derivative(3) = 1.0_8
            B_Functions_Derivative(4) = 0.0_8
            B_Functions_Derivative(5) = xi/(1.0_8 - zeta)
            B_Functions_Derivative(6) = 0.0_8
            B_Functions_Derivative(7) = 0.0_8
            B_Functions_Derivative(8) = (6.0_8*zeta - 1.0_8)
            B_Functions_Derivative(9) = xi*(6.0_8*zeta - 1.0_8)/(1.0_8 - zeta)
            B_Functions_Derivative(10) = 0.0_8
            B_Functions_Derivative(11) = 3.0_8*eta
            B_Functions_Derivative(12) = 0.5_8*(3.0_8*xi**2/(1.0_8 - zeta) - (1.0_8 - zeta))*(3.0_8*eta/(1.0_8 - zeta))
            B_Functions_Derivative(13) = 0.5_8*(3.0_8*xi**2 - (1.0_8 - zeta)**2)*1.0_8/(1.0_8 - zeta)
            B_Functions_Derivative(14) = (3.0_8*eta)*xi/(1.0_8 - zeta)

            Shape_Functions_Derivative(:) = matmul(VDM_Inverse, B_Functions_Derivative)

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Pyramidal_Shape_Functions_Derivative_eta

    Function Generate_Pyramidal_Shape_Functions_Derivative_zeta(this, Degree, a, eta, xi, zeta) Result(Value)

        type(PropertiesTypeD), intent(in) :: this

        real(kind = 8), dimension(:), allocatable :: Shape_Functions_Derivative
        real(kind = 8), dimension(:), allocatable :: B_Functions_Derivative
        integer, intent(in)                       :: a, Degree
        real(kind = 8)                            :: Value
        real(kind = 8), intent(in)                :: xi, eta, zeta
        real(kind = 8), dimension(:,:), allocatable :: VDM_Inverse

        Value = 0.0_8

        if(this%Isoparametric_Coordinates(1,1) == 10.0_8) print *, 'xi = 10'

        if (Degree == 1) then

            allocate(Shape_Functions_Derivative(5))

            Shape_Functions_Derivative(1) = 0.25_8*(-1.0_8 - (xi*eta)/(1.0_8 - zeta)**2)
            Shape_Functions_Derivative(4) = 0.25_8*(-1.0_8 + (xi*eta)/(1.0_8 - zeta)**2)
            Shape_Functions_Derivative(3) = 0.25_8*(-1.0_8 - (xi*eta)/(1.0_8 - zeta)**2)
            Shape_Functions_Derivative(2) = 0.25_8*(-1.0_8 + (xi*eta)/(1.0_8 - zeta)**2)
            Shape_Functions_Derivative(5) = 1.0_8

            Value = Shape_Functions_Derivative(a)

        else if (Degree == 2) then

            allocate(Shape_Functions_Derivative(14))
            allocate(B_Functions_Derivative(14))
            allocate(VDM_Inverse(14,14))

            VDM_Inverse(1,:) = (/ -0.020833333334583412_8, -0.027777777777777762_8, 0.027777777777777797_8, -0.030092592590509261_8, -0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555566_8, -0.055555555555555559_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(2,:) = (/ -0.020833333334583315_8, 0.027777777777777762_8, 0.027777777777777762_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, -0.055555555555555566_8, -0.055555555555555566_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, 0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(3,:) = (/ -0.020833333334583329_8, 0.027777777777777797_8, -0.027777777777777762_8, -0.030092592590509268_8, -0.16666666666666669_8, 0.018518518521851862_8, -0.055555555555555559_8, 0.055555555555555566_8, 0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, 0.16666666666666666_8 /)
            VDM_Inverse(4,:) = (/ -0.020833333334583315_8, -0.027777777777777797_8, -0.027777777777777797_8, -0.030092592590509261_8, 0.16666666666666669_8, 0.018518518521851858_8, 0.055555555555555559_8, 0.055555555555555559_8, -0.083333333333333329_8, 0.055555555555555552_8, 0.055555555555555552_8, 0.11111111111111110_8, -0.16666666666666666_8, -0.16666666666666666_8 /)
            VDM_Inverse(5,:) = (/ -0.050000000015000032_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.08333333358333331_8, 0.0000000000000000_8, 0.13333333337333336_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8, -0.0000000000000000_8 /)
            VDM_Inverse(6,:) = (/ 0.06666666666666652_8, 0.0_8, 0.27777777777777779_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.0_8, -0.055555555555555532_8, 0.0_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, -0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(7,:) = (/ 0.06666666666666652_8, 0.27777777777777779_8, 0.0000000000000000_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, -0.055555555555555539_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, -0.0000000000000000_8, -0.33333333333333331_8 /)
            VDM_Inverse(8,:) = (/ 0.06666666666666652_8, 0.0000000000000000_8, -0.27777777777777779_8, -0.037037037037037035_8, 0.0000000000000000_8, 0.0074074074074074112_8, 0.0000000000000000_8, 0.055555555555555539_8, -0.0000000000000000_8, -0.11111111111111110_8, 0.22222222222222221_8, -0.22222222222222221_8, 0.33333333333333331_8, -0.0000000000000000_8 /)
            VDM_Inverse(9,:) = (/ 0.06666666666666652_8, -0.27777777777777779_8, 0.0_8, -0.037037037037037035_8, 0.0_8, 0.0074074074074074112_8, 0.055555555555555532_8, 0.0_8, 0.0_8, 0.22222222222222221_8, -0.11111111111111110_8, -0.22222222222222221_8, 0.0000000000000000_8, 0.33333333333333331_8 /)
            VDM_Inverse(10,:) = (/ 0.15000000000500011_8, -0.16666666666666663_8, 0.16666666666666663_8, 0.08333333332500036_8, -0.16666666666666663_8, -0.06666666668000028_8, -0.16666666666666663_8, 0.16666666666666666_8, -0.16666666666666663_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(11,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(12,:) = (/ 0.15000000000499999_8, 0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, -0.16666666666666666_8, -0.06666666668000000_8, 0.16666666666666666_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(13,:) = (/ 0.15000000000499999_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.08333333332499999_8, 0.16666666666666666_8, -0.06666666668000000_8, -0.16666666666666666_8, -0.16666666666666666_8, 0.16666666666666666_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8 /)
            VDM_Inverse(14,:) = (/ 0.26666666666666661_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.14814814814814814_8, 0.0000000000000000_8, 0.029629629629629631_8, 0.0000000000000000_8, 0.0000000000000000_8, 0.0000000000000000_8, -0.44444444444444442_8, -0.44444444444444442_8, 0.44444444444444442_8, 0.0000000000000000_8, 0.0000000000000000_8 /)

            B_Functions_Derivative(1) = 0.0_8
            B_Functions_Derivative(2) = 0.0_8
            B_Functions_Derivative(3) = 0.0_8
            B_Functions_Derivative(4) = 4.0_8
            B_Functions_Derivative(5) = xi*eta/(1.0_8 - zeta)**2
            B_Functions_Derivative(6) = 30.0_8*zeta - 10.0_8
            B_Functions_Derivative(7) = xi*(6.0_8)
            B_Functions_Derivative(8) = eta*(6.0_8)
            B_Functions_Derivative(9) = xi*eta*(6.0_8)/(1.0_8 - zeta) + xi*eta*(6.0_8*zeta - 1.0_8)/(1.0_8 - zeta)**2
            B_Functions_Derivative(10) = (1.0_8 - zeta)
            B_Functions_Derivative(11) = (1.0_8 - zeta)
            B_Functions_Derivative(12) = 0.25_8*(3.0_8*xi**2/(1.0_8 - zeta)**2 + 1.0_8)*(3.0_8*eta**2/(1.0_8 - zeta) - (1.0_8 - zeta)) + 0.25_8*(3.0_8*xi**2/(1.0_8 - zeta) - (1.0_8 - zeta))*(3.0_8*eta**2/(1.0_8 - zeta)**2 + 1.0_8)
            B_Functions_Derivative(13) = eta + 0.5_8*(3.0_8*xi**2 - (1.0_8 - zeta)**2)*eta/(1.0_8 - zeta)**2
            B_Functions_Derivative(14) = 0.5_8*(3.0_8*xi*eta**2/(1.0_8 - zeta)**2 + xi)

            Shape_Functions_Derivative(:) = matmul(VDM_Inverse, B_Functions_Derivative)

            Value = Shape_Functions_Derivative(a)

        end if

    end Function Generate_Pyramidal_Shape_Functions_Derivative_zeta

    subroutine Calculate_Pyr_D_Matrix(Properties,N,i,D_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: j, k, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:,:) :: D_Matrix

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(dSFMatT(Num_Gauss_Points,Num_Nodes,3), dSFMat(Num_Gauss_Points,3,Num_Nodes))

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Pyramidal_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        D_Matrix = 0.0_8

        do k = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                dSFMat(j,1,k) = Generate_Pyramidal_Shape_Functions_Derivative_xi(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,2,k) = Generate_Pyramidal_Shape_Functions_Derivative_eta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))
                dSFMat(j,3,k) = Generate_Pyramidal_Shape_Functions_Derivative_zeta(Properties, N%Degree, k, eta(j), xi(j), zeta(j))

                dSFMatT(j,k,1) = dSFMat(j,1,k)
                dSFMatT(j,k,2) = dSFMat(j,2,k)
                dSFMatT(j,k,3) = dSFMat(j,3,k)

            end do

        end do

        do j = 1, Num_Gauss_Points

            D_Matrix = D_Matrix + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), (Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

        end do
        
    end subroutine Calculate_Pyr_D_Matrix

    subroutine Calculate_Pyr_Mass_Matrix(Properties,N,i,Mass_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, index_2, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), intent(inout) :: Mass_Matrix

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Pyramidal_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Mass_Matrix = 0.0_8

        do index_1 = 1, Num_Nodes

            do index_2 = 1, Num_Nodes

                do j = 1, Num_Gauss_Points

                    Mass_Matrix(index_1,index_2) = Mass_Matrix(index_1,index_2) + w(j)*Generate_Pyramidal_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*Generate_Pyramidal_Shape_Functions(Properties, N%Degree, index_2, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

                end do

            end do

        end do

    end subroutine Calculate_Pyr_Mass_Matrix

    subroutine Calculate_Pyr_Source_Vector(Properties,N,i,Source_Vector)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer             :: index_1, j, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), intent(inout) :: Source_Vector

        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w

        Num_Gauss_Points = Properties%Elements(i)%Number_of_Gauss_Points

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), zeta(Num_Gauss_Points), w(Num_Gauss_Points))

        call Generate_Pyramidal_Gauss_Points(Num_Gauss_Points, eta, xi, zeta, w)

        Source_Vector = 0.0_8

        do index_1 = 1, Num_Nodes

            do j = 1, Num_Gauss_Points

                Source_Vector(index_1) = Source_Vector(index_1) + w(j)*Generate_Pyramidal_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j), zeta(j))*ABS(Properties%Elements(i)%Det_Jacobian(j))

            end do

        end do

    end subroutine Calculate_Pyr_Source_Vector

    Function Pyr_Face_Boundary(Properties,N,i,Side_Nodes) Result(B_Matrix)

        type(PropertiesTypeD) :: Properties
        type(NTypeD)          :: N

        integer, intent(in) :: i
        integer, dimension(:) :: Side_Nodes
        integer             :: Num_Gauss_Points, a, b, index_1, index_2, k, l, p
        real(kind = 8), dimension(:,:), allocatable :: B_Matrix
        real(kind = 8), dimension(:), allocatable :: xi, eta, zeta, w
        integer, dimension(:), allocatable :: Nodes
        real(kind = 8), dimension(:,:), allocatable :: Shape_Functions
        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi, dz_dxi, dx_deta, dy_deta, dz_deta
        logical :: square

        allocate(B_Matrix(Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

        if (size(Side_Nodes) == (N%Degree+1)**2) then
            square = .true.
        else
            square = .false.
        end if

        if (square) then

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
                Nodes(5:9) = [6,7,8,9,14]
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
                    k = Side_Nodes(p)
                    Shape_Functions(k,l) = Generate_Pyramidal_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)
                    dx_dxi(l) = dx_dxi(l) + Generate_Pyramidal_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Pyramidal_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Pyramidal_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Pyramidal_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Pyramidal_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Pyramidal_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
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
                Nodes(4) = 5
                Nodes(5) = 6
                Nodes(6) = 7
            end if

            if (N%Degree == 2) then
                if (side_nodes(1) == 1 .and. side_nodes(2) == 2 .and. side_nodes(3) == 5) side_nodes = (/1,2,5,6,11,10/)
                if (side_nodes(1) == 2 .and. side_nodes(2) == 3 .and. side_nodes(3) == 5) side_nodes = (/2,3,5,7,12,11/)
                if (side_nodes(1) == 3 .and. side_nodes(2) == 4 .and. side_nodes(3) == 5) side_nodes = (/3,4,5,8,13,12/)
                if (side_nodes(1) == 1 .and. side_nodes(2) == 4 .and. side_nodes(3) == 5) side_nodes = (/1,4,5,9,13,10/)
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
                    k = Side_Nodes(p)
                    Shape_Functions(k,l) = Generate_Tetrahedral_Shape_Functions(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)
                    dx_dxi(l) = dx_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_dxi(l) = dy_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_dxi(l) = dz_dxi(l) + Generate_Tetrahedral_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
                    dx_deta(l) = dx_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,1)
                    dy_deta(l) = dy_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,2)
                    dz_deta(l) = dz_deta(l) + Generate_Tetrahedral_Shape_Functions_Derivative_eta(Properties, N%Degree, Nodes(p), eta(l), xi(l), 0.0_8)*Properties%Elements(i)%Coordinates(k,3)
                end do
    
            end do

        end if

        B_Matrix = 0.0_8

        do a = 1, size(Side_Nodes)

           do b = 1,  size(Side_Nodes)

                index_1 = Side_Nodes(a)
                index_2 = Side_Nodes(b)

                if (square) then
                    B_Matrix(index_1,index_2) = B_Matrix(index_1,index_2) + 0.5_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                else
                    if (index_1 /= 0 .and. index_2 /= 0) then
                        B_Matrix(index_1,index_2) = B_Matrix(index_1,index_2) + 0.25_8*sum(w*sqrt((dy_dxi*dz_deta - dz_dxi*dy_deta)**2 + (dz_dxi*dx_deta - dx_dxi*dz_deta)**2 + (dx_dxi*dy_deta - dy_dxi*dx_deta)**2)*Shape_Functions(index_1,:)*Shape_Functions(index_2,:))
                    end if
                end if

            end do

        end do

    end Function Pyr_Face_Boundary

end module