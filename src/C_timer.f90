Module C_timer
!
! PURPOSE:
!   To time functions so that the speed of the program can be calculated accuratly
!
! Record of Revsions:
!   Date        Programer       Change
!   04/11/22    J.I.M           Original code
!   23/06/23    J.I.M           Updated the timer code to use date
!
use m_Constants_Mod
implicit none
private
public t_timer

TYPE t_timer
    REAL(Kind=dp)   :: time
    REAL(Kind=dp)   :: start_time
    REAL(Kind=dp)   :: input_time
    REAL(Kind=dp)   :: matrix_time
    REAL(Kind=dp)   :: solve_time
    REAL(Kind=dp)   :: output_time
    integer         :: Secs,mins,milli,hours
CONTAINS
    PROCEDURE,pass :: Startdate
    PROCEDURE,pass :: elapseddate
    PROCEDURE,PASS :: start_timer
    PROCEDURE,PASS :: elapsed_time
    PROCEDURE,PASS :: input
    PROCEDURE,PASS :: matrix
    PROCEDURE,PASS :: solve
end type t_timer
    
contains

subroutine Startdate(this)
    class(t_timer)  :: this
    integer         :: values(8)
    call date_and_time(values=values)
    this%hours = values(5)
    this%mins = values(6)
    this%Secs = values(7)
    this%milli = values(8)
end subroutine

subroutine elapseddate(this)
    class(t_timer)  :: this
    integer         :: values(8)
    integer         :: newsec,newmin,newmilli,newhour
    call date_and_time(values=values)
    newhour = values(5)
    newmin = values(6)
    newsec = values(7)
    newmilli = values(8)

    newmilli = newmilli - this%milli
    if (newmilli < 0) then
        newmilli = 1000 + newmilli
        newsec = newsec - 1
    end if
    newsec = newsec - this%Secs
    if (newsec < 0) then
        newsec = 60 + newsec
        newmin = newmin - 1
    end if
    newmin = newmin - this%mins
    if (newmin < 0) then
        newmin = 60 + newmin
        newhour = newhour - 1
    end if
    newhour = newhour - this%hours
    if (newhour < 0) then
        newhour = 24 + newhour
        print*,"A DAY HAS PASSED"
    end if

    if (newmin < 10) then
        if (newsec < 10) then
            if (newmilli < 100) then
                if (newmilli < 10) then
            write(*,'(a,i2,a,a,i1,a,a,i1,a,i1)')" THE ELASED TIME HAS BEEN:", newhour,"::","0",newmin,":","0",newsec,".00",newmilli
                else
            write(*,'(a,i2,a,a,i1,a,a,i1,a,i2)') " THE ELASED TIME HAS BEEN:",newhour,"::","0",newmin,":","0",newsec,".0",newmilli
                end if
            else
            write(*,'(a,i2,a,a,i1,a,a,i1,a,i3)')" THE ELASED TIME HAS BEEN:", newhour,"::","0",newmin,":","0",newsec,".",newmilli
            end if
        else
            if (newmilli < 100) then
                if (newmilli < 10) then
                write(*,'(a,i2,a,a,i1,a,i2,a,i1)') " THE ELASED TIME HAS BEEN:",newhour,"::","0",newmin,":",newsec,".00",newmilli
                else
                    write(*,'(a,i2,a,a,i1,a,i2,a,i2)')" THE ELASED TIME HAS BEEN:", newhour,"::","0",newmin,":",newsec,".0",newmilli
                end if
            else
                write(*,'(a,i2,a,a,i1,a,i2,a,i3)') " THE ELASED TIME HAS BEEN:",newhour,"::","0",newmin,":",newsec,".",newmilli
            end if
        end if
    else
        if (newsec < 10) then
            if (newmilli < 100) then
                if (newmilli < 10) then
                write(*,'(a,i2,a,i2,a,a,i1,a,i1)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":","0",newsec,".00",newmilli
                else
                    write(*,'(a,i2,a,i2,a,a,i1,a,i2)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":","0",newsec,".0",newmilli
                end if
            else
                write(*,'(a,i2,a,i2,a,a,i1,a,i3)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":","0",newsec,".",newmilli
            end if
        else
            if (newmilli < 100) then
                if (newmilli < 10) then
                    write(*,'(a,i2,a,i2,a,i2,a,i1)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":",newsec,".00",newmilli
                else
                    write(*,'(a,i2,a,i2,a,i2,a,i2)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":",newsec,".0",newmilli
                end if
            else
                write(*,'(a,i2,a,i2,a,i2,a,i3)') " THE ELASED TIME HAS BEEN:",newhour,"::",newmin,":",newsec,".",newmilli
            end if
        end if
    end if
    
    this%hours = values(5)
    this%mins = values(6)
    this%Secs = values(7)
    this%milli = values(8)
end subroutine

subroutine start_timer(this) !starts the timer
    CLASS(t_timer),INTENT(INOUT) :: this
    CALL CPU_TIME(this%time)
    this%start_time = this%time
end subroutine start_timer

subroutine input(this) !calls the time taken
    CLASS(t_timer),INTENT(INOUT) :: this
    real(kind=dp) :: time
    CALL CPU_TIME(time)
    this%input_time = time - this%time
    this%time = time
end subroutine input

subroutine matrix(this) !calls the time taken
    CLASS(t_timer),INTENT(INOUT) :: this
    real(kind=dp) :: time
    CALL CPU_TIME(time)
    this%matrix_time = time - this%time
    this%time = time
end subroutine matrix

subroutine solve(this) !calls the time taken
    CLASS(t_timer),INTENT(INOUT) :: this
    real(kind=dp) :: time
    CALL CPU_TIME(time)
    this%solve_time = time - this%time
    this%time = time
end subroutine solve

subroutine elapsed_time(this) !calls the time taken
    CLASS(t_timer),INTENT(INOUT) :: this
    REAL(KIND = dp) :: stopTime
    CALL CPU_TIME(stopTime)

    this%output_time = stopTime - this%time

    print *, 'Stage           ', 'Total CPU Time', '            Percent CPU Time'
    print *, '------------------------------------------------------------'
    print *, "1   Total Run", (stopTime - this%start_time), 100.0
    print *, "2   Input    ", this%input_time, (this%input_time/(stopTime - this%start_time))*100
    print *, "3   Matrix   ", this%matrix_time, (this%matrix_time/(stopTime - this%start_time))*100
    print *, "4   Solve    ", this%solve_time, (this%solve_time/(stopTime - this%start_time))*100
    print *, "5   Output   ", this%output_time, (this%output_time/(stopTime - this%start_time))*100
    print *, " "
    
end subroutine elapsed_time

end module