module Determine_Solve_Method
!
! Purpose:
! To read from an input file and determine whether the problem is a transport problem or a diffusion problem.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
implicit none
    
contains

    subroutine Determine_Solve(Solve_Transport_Method)

        logical :: Solve_Transport_Method

        integer :: Angle

        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) Angle

        if (Angle == 0) then

            Solve_Transport_Method = .false.

        else if (Angle > 0) then

            Solve_Transport_Method = .true.

        else

            write(*,*) 'Invalid Angle'
            stop

        end if

        rewind(1)

    end subroutine Determine_Solve
    
end module Determine_Solve_Method