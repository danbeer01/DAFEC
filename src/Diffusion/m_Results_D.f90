module m_Results_D
!
! Purpose:
! To define a derived type for the results of the simulation and a procedure to print the k_eff
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
implicit none

type :: ResultsTypeD
    real(kind=8), dimension(:,:), allocatable  :: Flux
    real(kind=8)                               :: k_eff
CONTAINS
    PROCEDURE, pass :: print_keff_Diffusion
end type ResultsTypeD

contains

    subroutine print_keff_Diffusion(Results)
        class(ResultsTypeD)      :: Results
        
        print *, "k_eff = ", Results%k_eff
        print *, " "
    end subroutine print_keff_Diffusion

end module