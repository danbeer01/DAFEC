module m_Results
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

type :: ResultsType
    real(kind=8), dimension(:,:,:), allocatable  :: Angular_Flux
    real(kind=8), dimension(:,:), allocatable    :: Scalar_Flux
    real(kind=8)                                 :: k_eff
CONTAINS
    PROCEDURE, pass :: print_keff
end type ResultsType

contains

    subroutine print_keff(Results)
        class(ResultsType)      :: Results

        write(*, '(A)') " **************************************"

        write(*, '(A)') " *                                    *"

        write(*, *)      "* ", "K-eff = ", Results%k_eff,   "*"

        write(*, '(A)') " *                                    *"

        write(*, '(A)') " **************************************"

        print *, " "

    end subroutine print_keff

end module m_Results