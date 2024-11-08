module m_Sweep_Order_1D
!
! Purpose:
! To determine the order in which elements are swept in a 1D transport problem.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
    use m_Read_Properties
    use m_Calculate_mu_w

    implicit none

    contains

    subroutine determine_sweep_order_1D(Properties, N, Sweep_Order)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in) :: N

        integer, intent(out), dimension(:,:) :: Sweep_Order

        real(8), dimension(N%Ordinates) :: mu

        integer :: i, ang, side_index

        call calculate_mu_1D(mu)

        do ang = 1, N%Angle

            if (mu(ang) > 0.0) then

                do i = 1, N%Element

                    Sweep_Order(ang,i) = i

                end do

            else

                do i = 1, N%Element

                    Sweep_Order(ang,i) = N%Element + 1 - i

                end do

            end if

        end do

        do i = 1, N%Element

            Properties%Elements(i)%Neighbours(1,1) = i - 1

            Properties%Elements(i)%Neighbours(2,1) = i + 1

            if (i == N%Element) Properties%Elements(i)%Neighbours(2,1) = 0

            do side_index = 1,2

                allocate(Properties%Elements(i)%Sides(side_index)%F_in_Matrix(N%Angle,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))

                allocate(Properties%Elements(i)%Sides(side_index)%Boundary(N%Group,N%Ordinates,Properties%Elements(i)%Number_of_Nodes))

            end do

        end do

    end subroutine determine_sweep_order_1D

end module m_Sweep_Order_1D

