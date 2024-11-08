module m_Boundary_Conditions
!
! Purpose:
! To calculate the boundary conditions for the problem.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties

implicit none

contains

    subroutine Boundary_Conditions_1D(Properties, N, i, j, k, ang, Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: i, j, k, ang

        real(kind = 8), dimension(:), intent(inout) :: Boundary

        Boundary = 0.0_8

        if (i == 1 .and. j == 1) then

            if (Properties%LBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,N%Angle+1-ang,:)

            else if (Properties%LBC == 2) then

                Boundary = 0.0_8

            else if (Properties%LBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,N%Angle+1-ang,:)

            else if (Properties%LBC == 4) then

                Boundary = Properties%Elements(i)%Flux(k,ang,:)

            else if (Properties%LBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (i == N%Element .and. j == 2) then

            if (Properties%RBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,N%Angle+1-ang,:)

            else if (Properties%RBC == 2) then

                Boundary = 0.0_8

            else if (Properties%RBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,N%Angle+1-ang,:)

            else if (Properties%RBC == 4) then

                Boundary = Properties%Elements(i)%Flux(k,ang,:)

            else if (Properties%RBC == 5) then

                Boundary = Properties%Q_s

            end if

        end if

    end subroutine Boundary_Conditions_1D

    subroutine Boundary_Conditions_2D(Properties, N, i, j, k, ang, Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: i, j, k, ang

        integer :: reflect_counter_1, reflect_counter_2

        real(kind = 8), dimension(:), intent(inout) :: Boundary

        Boundary = 0.0_8

        reflect_counter_1 = 0
        reflect_counter_2 = 0

        ! Counters for reflective boundary conditions

        if (ang <= N%Ordinates/4) then

            reflect_counter_1 = ang + N%Ordinates/4

            reflect_counter_2 = ang + N%Ordinates/2

        else if (ang <= N%Ordinates/2) then

            reflect_counter_1 = ang - N%Ordinates/4

            reflect_counter_2 = ang + N%Ordinates/2

        else if (ang <= 3*N%Ordinates/4) then

            reflect_counter_1 = ang + N%Ordinates/4

            reflect_counter_2 = ang - N%Ordinates/2

        else if (ang <= N%Ordinates) then

            reflect_counter_1 = ang - N%Ordinates/4

            reflect_counter_2 = ang - N%Ordinates/2

        end if

        ! Now implement the boundary conditions

        if (Properties%Elements(i)%Neighbours(j,2) == 1) then

            if (Properties%LBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%LBC == 2) then

                Boundary = 0.0_8

            else if (Properties%LBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%LBC == 4) then

                print *, 'Periodic Boundaries not implemented'

            else if (Properties%LBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 2) then

            if (Properties%RBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%RBC == 2) then

                Boundary = 0.0_8

            else if (Properties%RBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%RBC == 4) then

                print *, 'Periodic Boundaries not implemented'

            else if (Properties%RBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 3) then

            if (Properties%TBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%TBC == 2) then

                Boundary = 0.0_8

            else if (Properties%TBC == 3) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%TBC == 4) then

                print *, 'Periodic Boundaries not implemented'

            else if (Properties%TBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 4) then

            if (Properties%BBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%BBC == 2) then

                Boundary = 0.0_8

            else if (Properties%BBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%BBC == 4) then

                print *, 'Periodic Boundaries not implemented'

            else if (Properties%BBC == 5) then

                Boundary = Properties%Q_s

            end if

        end if

    end subroutine Boundary_Conditions_2D

    subroutine Boundary_Conditions_3D(Properties, N, i, j, k, ang, Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: i, j, k, ang

        integer :: reflect_counter_1, reflect_counter_2, reflect_counter_3

        real(kind = 8), dimension(:), intent(inout) :: Boundary

        Boundary = 0.0_8

        reflect_counter_1 = 0
        reflect_counter_2 = 0
        reflect_counter_3 = 0

        ! Counters for reflective boundary conditions

        if (ang <= N%Ordinates/8) then

            reflect_counter_1 = ang + N%Ordinates/8

            reflect_counter_2 = ang + N%Ordinates/4

            reflect_counter_3 = ang + N%Ordinates/2

        else if (ang <= N%Ordinates/4) then

            reflect_counter_1 = ang - N%Ordinates/8

            reflect_counter_2 = ang + N%Ordinates/4

            reflect_counter_3 = ang + N%Ordinates/2

        else if (ang <= 3*N%Ordinates/8) then

            reflect_counter_1 = ang + N%Ordinates/8

            reflect_counter_2 = ang - N%Ordinates/4

            reflect_counter_3 = ang + N%Ordinates/2

        else if (ang <= N%Ordinates/2) then

            reflect_counter_1 = ang - N%Ordinates/8

            reflect_counter_2 = ang - N%Ordinates/4

            reflect_counter_3 = ang + N%Ordinates/2

        else if (ang <= 5*N%Ordinates/8) then

            reflect_counter_1 = ang + N%Ordinates/8

            reflect_counter_2 = ang + N%Ordinates/4

            reflect_counter_3 = ang - N%Ordinates/2

        else if (ang <= 3*N%Ordinates/4) then

            reflect_counter_1 = ang - N%Ordinates/8

            reflect_counter_2 = ang + N%Ordinates/4

            reflect_counter_3 = ang - N%Ordinates/2

        else if (ang <= 7*N%Ordinates/8) then

            reflect_counter_1 = ang + N%Ordinates/8

            reflect_counter_2 = ang - N%Ordinates/4

            reflect_counter_3 = ang - N%Ordinates/2

        else if (ang <= N%Ordinates) then

            reflect_counter_1 = ang - N%Ordinates/8

            reflect_counter_2 = ang - N%Ordinates/4

            reflect_counter_3 = ang - N%Ordinates/2

        end if

        ! Now implement the boundary conditions

        if (Properties%Elements(i)%Neighbours(j,2) == 1) then

            if (Properties%LBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%LBC == 2) then

                Boundary = 0.0_8

            else if (Properties%LBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%LBC == 4) then

                ! Boundary = Properties%Elements(1)%Flux(k,ang,1)

            else if (Properties%LBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 2) then

            if (Properties%RBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%RBC == 2) then

                Boundary = 0.0_8

            else if (Properties%RBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_1,:)

            else if (Properties%RBC == 4) then

                ! Boundary = Properties%Elements(N%Element)%Flux(k,ang,size(Properties%Elements(N%Element)%Flux,3))

            else if (Properties%RBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 3) then

            if (Properties%TBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%TBC == 2) then

                Boundary = 0.0_8

            else if (Properties%TBC == 3) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%TBC == 4) then

                ! Boundary = Properties%Elements(1)%Flux(k,N%Angle+1-ang,1)

            else if (Properties%TBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 4) then

            if (Properties%BBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%BBC == 2) then

                Boundary = 0.0_8

            else if (Properties%BBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_2,:)

            else if (Properties%BBC == 4) then

                ! Boundary = Properties%Elements(N%Element)%Flux(k,N%Angle+1-ang,size(Properties%Elements(N%Element)%Flux,3))

            else if (Properties%BBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 5) then

            if (Properties%FBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_3,:)

            else if (Properties%FBC == 2) then

                Boundary = 0.0_8

            else if (Properties%FBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_3,:)

            else if (Properties%FBC == 4) then

                ! Boundary = Properties%Elements(N%Element)%Flux(k,ang,size(Properties%Elements(N%Element)%Flux,3))

            else if (Properties%FBC == 5) then

                Boundary = Properties%Q_s

            end if

        else if (Properties%Elements(i)%Neighbours(j,2) == 6) then

            if (Properties%BaBC == 1) then

                Boundary = Properties%Elements(i)%Flux(k,reflect_counter_3,:)

            else if (Properties%BaBC == 2) then

                Boundary = 0.0_8

            else if (Properties%BaBC == 3) then

                Boundary = Properties%alpha*Properties%Elements(i)%Flux(k,reflect_counter_3,:)

            else if (Properties%BaBC == 4) then

                ! Boundary = Properties%Elements(N%Element)%Flux(k,ang,size(Properties%Elements(N%Element)%Flux,3))

            else if (Properties%BaBC == 5) then

                Boundary = Properties%Q_s

            end if

        end if

    end subroutine Boundary_Conditions_3D

end module m_Boundary_Conditions