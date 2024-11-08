module m_Solver_D
!
! Purpose:
! To solve the neutron diffusion equation using the finite element method
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Results_D
use m_PETSc
use m_Construct_Matrix_D
use m_VTK_Reader
use m_Normalisation_D

implicit none

contains

    subroutine Solver_Diffusion(Properties, Results, N, this)

        type(PropertiesTypeD), intent(inout) :: Properties
        type(ResultsTypeD), intent(inout)    :: Results
        type(NTypeD), intent(in)             :: N
        type(PETScGroupType)                 :: PETSc
        type(Mesh), intent(inout)            :: this

        real(kind = 8), dimension(N%Group,N%N)      :: Total_Source, Total_Source_new
        real(kind = 8), dimension(N%Group,N%N)      :: Flux
        real(kind = 8)                              :: lambda_old = 1.1_8, lambda_new = 1.0_8

        integer                                     :: max_iter = 100000 ! maximum number of iterations
        real(kind = 8)                              :: tol = 1.0e-8 ! tolerance for convergence

        integer                                     :: iter = 0

        integer, dimension(:,:), allocatable        :: Periodic_Pairs

        integer :: i, g_index

        write(*, '(A)') "*******************************"

        write(*, '(A)') "*                             *"

        write(*, '(A, A, A)') "*       ", "DIFFUSION SOLVE ", "      *"

        write(*, '(A)') "*                             *"

        write(*, '(A)') "*******************************"
        
        ! Allocate memory for the flux array
        allocate(Results%Flux(N%Group,N%N))

        ! Set the initial flux guesses and initialise the sources
        do i = 1,N%Element

            Flux = 1.0_8

            Properties%Elements(i)%Source = 0.0_8

        end do

        ! Construct the matrix
        call Construct_Matrix(Properties, N, Periodic_Pairs)

        ! Initialise PETSc
        call Start_PETSc(PETSc, Properties, this, N, tol, max_iter, Periodic_Pairs)

        ! Iterate until the eigenvalue converges
        do while (ABS((lambda_new - lambda_old)/(lambda_old)) > tol)

            lambda_old = lambda_new

            call Calculate_Element_Flux(Properties, N, Flux)

            call Calculate_Source(Properties, N, lambda_new, Total_Source)

            do g_index = 1,N%Group

                call KSP_Solve(PETSc, Total_Source(g_index,:), N%N, Flux(g_index,:), g_index)

            end do

            call Calculate_Element_Flux(Properties, N, Flux)

            call Calculate_Source(Properties, N, lambda_new, Total_Source_new)

            lambda_new=lambda_old*((SUM(Total_Source_new))/(SUM(Total_Source))) ! Calculate the new eigenvalue

            iter = iter + 1

            print *, 'Iteration ', iter, 'Eigenvalue', lambda_new, 'Error: ', ABS((lambda_new - lambda_old)/(lambda_old))

        end do

        print *, 'Iterations: ', iter

        ! End PETSc
        call End_PETSc(PETSc)

        Results%Flux = Flux

        if (Properties%Case == 1) then

            Results%k_eff = lambda_new

            ! Normalise the flux
            call Normalisation(Properties, Results, N)

        end if

    end subroutine Solver_Diffusion

    subroutine Calculate_Source(Properties, N, k_eff, Total_Source)

        type(PropertiesTypeD), intent(inout) :: Properties
        type(NTypeD), intent(in)          :: N

        real(kind = 8), dimension(:,:) :: Total_Source

        real(kind = 8) :: k_eff

        integer :: i, a, g_index, g_s_index

        do i = 1,N%Element

            do g_index = 1,N%Group

                if (Properties%Case == 0) then

                    Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Sigma_f(g_index)*Properties%Elements(i)%Source_Vector(:)

                else if (Properties%Case == 1) then

                    Properties%Elements(i)%Source(g_index,:) = (Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_index,:))

                end if

                do g_s_index = g_index,2,-1 ! Upscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index-1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index-1,:))

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index-1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index-1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index-1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index-1,g_index) + Properties%Chi(g_s_index-1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index-1,:))

                    end if

                end do

                do g_s_index = g_index,(N%Group-1) ! Downscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index+1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index+1,:))

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index+1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index+1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index+1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,:) = Properties%Elements(i)%Source(g_index,:) + (Properties%Elements(i)%Sigma_s(g_s_index+1,g_index) + Properties%Chi(g_s_index+1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(g_s_index+1,:))

                    end if

                end do

            end do

        end do

        Total_Source = 0.0_8

        do i = 1, N%Element

            do g_index = 1,N%Group

                do a = 1, Properties%Elements(i)%Number_of_Nodes

                    Total_Source(g_index,Properties%Elements(i)%Cell_Pointers(a)) = Total_Source(g_index,Properties%Elements(i)%Cell_Pointers(a)) + Properties%Elements(i)%Source(g_index,a)

                end do

            end do

        end do

    end subroutine Calculate_Source

    subroutine Calculate_Element_Flux(Properties, N, Flux)

        type(PropertiesTypeD), intent(inout) :: Properties
        type(NTypeD), intent(in)          :: N

        real(kind = 8), dimension(:,:) :: Flux

        integer :: i, a, g_index

        do i = 1,N%Element

            do g_index = 1,N%Group

                do a = 1, Properties%Elements(i)%Number_of_Nodes

                    Properties%Elements(i)%Flux(g_index,a) = Flux(g_index,Properties%Elements(i)%Cell_Pointers(a))

                end do

            end do

        end do

    end subroutine Calculate_Element_Flux

end module