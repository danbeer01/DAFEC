module m_normalisation_D
!
! Purpose:
! To normalise the flux solution
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Constants_Mod
use m_Read_Properties_D
use m_Results_D
use m_VTK_Reader
use m_Gauss_Points

implicit none

contains

    subroutine normalisation(Properties, Results, N)

    type(PropertiesTypeD), intent(inout) :: Properties
    type(ResultsTypeD), intent(inout)    :: Results
    type(NTypeD), intent(in)             :: N

    real(kind=8) :: Norm, R

    integer :: k, i

    Norm = 0

    do i = 1,N%Element

        if (Properties%g == 1) then

            R = SUM(Properties%Elements(i)%Coordinates(:,1))/Properties%Elements(i)%Number_of_Nodes

            Properties%Elements(i)%Volume = 2.0_8*PI*R*Properties%Elements(i)%Volume

        else if (Properties%g == 2) then

            R = SUM(Properties%Elements(i)%Coordinates(:,1))/Properties%Elements(i)%Number_of_Nodes

            Properties%Elements(i)%Volume = 4.0_8*PI*(R**2)*Properties%Elements(i)%Volume

        end if

        do k = 1,N%Group

            if (.not. Properties%Adjoint) then

                Norm = Norm + (1.0_8/Results%k_eff)*Properties%Elements(i)%Sigma_f(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(k,:)))

            else if (Properties%Adjoint) then

                Norm = Norm + sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(k,:)))

            end if

        end do

    end do

    if (Properties%g == 1) Norm = Norm*2.0_8*PI
    if (Properties%g == 2) Norm = Norm*4.0_8*PI

    Results%Flux = Results%Flux / Norm

    do i = 1, N%Element

        Properties%Elements(i)%Flux = Properties%Elements(i)%Flux / Norm

    end do
    
    end subroutine normalisation

end module