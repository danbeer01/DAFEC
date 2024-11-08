module m_normalisation
!
! Purpose:
! To normalise the Scalar_Flux solution
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Constants_Mod
use m_Read_Properties
use m_Results
use m_VTK_Reader
use m_Gauss_Points

implicit none

contains

    subroutine normalisation(Properties, Results, N)

    type(PropertiesType), intent(inout) :: Properties
    type(ResultsType), intent(inout) :: Results
    type(NType), intent(in)          :: N

    real(kind=8) :: Norm

    real(kind=8) :: R

    real(kind=8), dimension(N%Group,N%Element) :: Element_Scalar_Flux

    integer :: k, i

    Norm = 0

    Element_Scalar_Flux = 0.0_8

    do k = 1,N%Group

        do i = 1,N%Element

            Element_Scalar_Flux(k,i) = sum(Properties%Elements(i)%Scalar_Flux(k,:))

            Element_Scalar_Flux(k,i) = Element_Scalar_Flux(k,i) / Properties%Elements(i)%Number_of_Nodes

        end do

    end do
    
    do i = 1,N%Element

        if (Properties%g == 1) then

            R = SUM(Properties%Elements(i)%Coordinates(:,1))/Properties%Elements(i)%Number_of_Nodes

            Properties%Elements(i)%Volume = 2.0_8*PI*R*Properties%Elements(i)%Volume

        else if (Properties%g == 2) then

            R = SUM(Properties%Elements(i)%Coordinates(:,1))/Properties%Elements(i)%Number_of_Nodes

            Properties%Elements(i)%Volume = 4.0_8*PI*(R**2)*Properties%Elements(i)%Volume

        end if

        do k = 1,N%Group

            if (Properties%Adjoint == 0) then

                Norm = Norm + (1.0_8/Results%k_eff)*Properties%Elements(i)%Sigma_f(k)*Element_Scalar_Flux(k,i)*Properties%Elements(i)%Volume

            else if (Properties%Adjoint == 1) then

                Norm = Norm + Element_Scalar_Flux(k,i)*Properties%Elements(i)%Volume

            end if

        end do

    end do

    Results%Scalar_Flux = Results%Scalar_Flux / Norm

    Results%Angular_Flux = Results%Angular_Flux / Norm

    do i = 1,N%Element

        Properties%Elements(i)%Scalar_Flux = Properties%Elements(i)%Scalar_Flux / Norm

        Properties%Elements(i)%Flux = Properties%Elements(i)%Flux / Norm

    end do
    
    end subroutine normalisation

end module m_normalisation