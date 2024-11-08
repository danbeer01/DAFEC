module m_Solver
!
! Purpose:
! To solve the neutron diffusion equation using the finite element method
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties
use m_Results
use m_Construct_Matrix_1D
use m_Construct_Matrix_2D
use m_Construct_Matrix_3D
use m_Create_Hexahedral_Shape_Functions
use m_Create_Tetrahedral_Shape_Functions
use m_Create_Prismatic_Shape_Functions
use m_Create_Pyramidal_Shape_Functions
use m_Create_Quadrilateral_Shape_Functions
use m_Create_Triangular_Shape_Functions
use m_Create_Shape_Functions
use m_VTK_Reader
use m_Normalisation
use m_Calculate_mu_w
use m_Small_Matrix_Solver
use m_Sweep_Order_1D
use m_Sweep_Order_2D
use m_Sweep_Order_3D
use m_Boundary_Conditions
use m_Cylindrical
use m_Spherical
use m_RZ

implicit none

contains

    subroutine Solver(Properties, Results, N)

        type(PropertiesType), intent(inout) :: Properties
        type(ResultsType), intent(inout)    :: Results
        type(NType), intent(in)             :: N

        real(kind = 8)                       :: Total_Source, Total_Source_new
        real(kind = 8)                       :: lambda_old = 1.1_8, lambda_new = 1.0_8

        integer, dimension(N%Ordinates,N%Element) :: Sweep_Order

        integer                             :: max_iter = 10000 ! maximum number of iterations
        real(kind = 8)                      :: tol = 1.0e-8 ! tolerance for convergence

        integer :: i, g_index, ang, iter, element_index

        logical :: solve = .true.

        write(*, '(A)') "*******************************"

        write(*, '(A)') "*                             *"

        write(*, '(A, A, A)') "*       ", "TRANSPORT SOLVE ", "      *"

        write(*, '(A)') "*                             *"

        write(*, '(A)') "*******************************"

        ! Allocate memory for the flux array
        allocate(Results%Angular_Flux(N%Group,N%N,N%Ordinates))
        allocate(Results%Scalar_Flux(N%Group,N%N))

        if (N%D == 1) call Calculate_Isoparametric_Coordinates(N%Degree, Properties)
        if (N%D == 2) call Calculate_Isoparametric_Quadrilateral_Coordinates(N%Degree, Properties)
        if (N%D == 3) call Calculate_Isoparametric_Hexahedral_Coordinates(N%Degree, Properties)

        if (N%D == 1) call Determine_Sweep_Order_1D(Properties, N, Sweep_Order)
        if (N%D == 2) call Determine_Sweep_Order_2D(Properties, N, Sweep_Order)
        if (N%D == 3) call Determine_Sweep_Order_3D(Properties, N, Sweep_Order)

        ! Set the initial flux guesses and initialise the sources
        do i = 1,N%Element

            Properties%Elements(i)%Scalar_Flux = 1.0_8

            Properties%Elements(i)%Flux = 1.0_8

            Properties%Elements(i)%Legendre_Flux = 0.0_8

            Properties%Elements(i)%Total_Source = 0.0_8

            Properties%Elements(i)%Sweep_Source = 0.0_8

            Properties%Elements(i)%Source = 0.0_8

        end do

        call Create_Matrices(Properties, N)

        iter = 0

        if (N%D == 1) then

            if (Properties%g == 1) then

                call Cylindrical_Solver(Properties, N, lambda_new)
                solve = .false.

            else if (Properties%g == 2) then

                call Spherical_Solver(Properties, N, lambda_new)
                solve = .false.

            end if

        else if (N%D == 2) then

            if (Properties%g == 1) then

                call RZ_Solver(Properties, N, Sweep_Order, lambda_new)
                solve = .false.

            end if

        end if

        if(solve) then

        ! Iterate until the eigenvalue converges
        do while (ABS((lambda_new - lambda_old)/(lambda_old)) > tol)

            lambda_old = lambda_new

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Source(Properties, N, lambda_new, i)

                if(N%Anisotropy > 0) call Anisotropic_Scattering(Properties, N, i)

            end do

            call Calculate_Total_Source(Total_Source, N, Properties)

            call Calculate_Boundaries(Properties, N)

            !$OMP PARALLEL DO PRIVATE(ang,i,element_index)
            do ang = 1,N%Ordinates

                do g_index = 1,N%Group

                    do i = 1,N%Element

                        element_index = Sweep_Order(ang,i)

                        call Up_Wind_Source(Properties, g_index, ang, element_index)

                        call Solve_Matrix(Properties%Elements(element_index)%K_Matrix(g_index,ang,:,:), Properties%Elements(element_index)%Total_Source(g_index,ang,:), Properties%Elements(element_index)%Flux(g_index,ang,:))

                    end do

                end do

            end do
            !$OMP END PARALLEL DO

            call Calculate_Scalar_Flux(Properties, N)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Source(Properties, N, lambda_new, i)

                if(N%Anisotropy > 0) call Anisotropic_Scattering(Properties, N, i)

            end do

            call Calculate_Total_Source(Total_Source_new, N, Properties)

            lambda_new=lambda_old*Total_Source_new/Total_Source ! Calculate the new eigenvalue

            iter = iter + 1

            if (iter > max_iter) then

                write(*,*) 'Maximum number of iterations reached'
                exit

            end if

            print *, 'Iteration ', iter, 'Eigenvalue', lambda_new, 'Error: ', ABS((lambda_new - lambda_old)/(lambda_old))

        end do

        print *, 'Iterations', iter

        end if

        call Input_Results(Results, Properties, N)

        if (Properties%Case == 1) then

            Results%k_eff = lambda_new

            ! Normalise the flux
            call Normalisation(Properties, Results, N)

        end if

    end subroutine Solver

    subroutine Calculate_Source(Properties, N, k_eff, i)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        real(kind = 8) :: k_eff

        integer, intent(in) :: i
        integer :: g_s_index, g_index, ang

        do ang = 1,N%Ordinates

            do g_index = 1,N%Group

                if (Properties%Case == 0) then

                    Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Sigma_f(g_index)*Properties%Elements(i)%Source_Vector(:) + Properties%Elements(i)%Sigma_s(0,g_index,g_index)*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_index,:))

                else if (Properties%Case == 1) then

                    Properties%Elements(i)%Source(g_index,ang,:) = (Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index) + Properties%Elements(i)%Sigma_s(0,g_index,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_index,:))

                end if

                do g_s_index = g_index,2,-1 ! Upscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index-1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_s_index-1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                    end if

                end do

                do g_s_index = g_index,(N%Group-1) ! Downscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index+1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_s_index+1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                    end if

                end do

            end do

        end do

    end subroutine Calculate_Source

    subroutine Anisotropic_Scattering(Properties, N, i)

        type(NType), intent(in)             :: N
        type(PropertiesType), intent(inout) :: Properties

        real(kind = 8) :: Y_l_m, P_l, mu_val, eta_val, xsi_val

        real(kind = 8), dimension(N%Angle/2) :: mu
        real(kind = 8), dimension(N%Angle)   :: mu_1D
        real(kind = 8), dimension(N%Ordinates/8) :: w

        integer, intent(in) :: i
        integer :: ang, g_s_index, g_index, w_index
        integer :: l, m

        if (N%D == 1) then

        call calculate_mu_1D(mu_1D)

        call calculate_w_1D(w)

        do g_index = 1, N%Group

            do l = 1, N%Anisotropy

                do ang = 1, N%Angle

                    call calculate_P_l(mu_1D(ang), l, P_l)

                    Properties%Elements(i)%Legendre_Flux(g_index,l,1,:) = Properties%Elements(i)%Legendre_Flux(g_index,l,1,:) + 0.5_8*w(ang)*P_l*Properties%Elements(i)%Flux(g_index,ang,:)

                end do

            end do

        end do

        do g_index = 1, N%Group

            do ang = 1, N%Angle

                do l = 1, N%Anisotropy

                    call calculate_P_l(mu_1D(ang), l, P_l)

                    Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (2*l+1)*P_l*Properties%Elements(i)%Sigma_s(l,g_index,g_index)*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_index,l,1,:))

                    do g_s_index = g_index,2,-1 ! Upscatter

                        Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (2*l+1)*P_l*(Properties%Elements(i)%Sigma_s(l,g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_s_index-1,l,1,:))

                    end do

                    do g_s_index = g_index,(N%Group-1) ! Downscatter

                        Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (2*l+1)*P_l*(Properties%Elements(i)%Sigma_s(l,g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_s_index+1,l,1,:))

                    end do

                end do

            end do

        end do

        else

        call calculate_mu(mu)

        call calculate_w(w)

        do g_index = 1, N%Group

            do l = 1, N%Anisotropy

                do m = 0, l

                    do ang = 1, N%Ordinates

                        if (N%D == 2) w_index = MOD(ang,N%Ordinates/4)
                        if (N%D == 3) w_index = MOD(ang,N%Ordinates/8)

                        if (w_index == 0) then

                            if (N%D == 2) w_index = N%Ordinates/4
                            if (N%D == 3) w_index = N%Ordinates/8

                        end if

                        call Ordinates(ang, N%Angle, mu, mu_val, eta_val, xsi_val)

                        Y_l_m = Calculate_Y_l_m(mu_val,eta_val,xsi_val,l,m)

                        Properties%Elements(i)%Legendre_Flux(g_index,l,m,:) = Properties%Elements(i)%Legendre_Flux(g_index,l,m,:) + (1.0_8/(2.0_8**N%D))*w(w_index)*Y_l_m*Properties%Elements(i)%Flux(g_index,ang,:)

                    end do

                end do

            end do

        end do

        do g_index = 1, N%Group

            do ang = 1, N%Ordinates

                call Ordinates(ang, N%Angle, mu, mu_val, eta_val, xsi_val)

                do l = 1, N%Anisotropy

                    do m = 0, l

                        Y_l_m = Calculate_Y_l_m(mu_val,eta_val,xsi_val,l,m)

                        Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + REAL((2.0_8*l + 1),8)*Y_l_m*Properties%Elements(i)%Sigma_s(l,g_index,g_index)*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_index,l,m,:))

                        do g_s_index = g_index,2,-1 ! Upscatter

                            Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + REAL((2.0_8*l + 1),8)*Y_l_m*(Properties%Elements(i)%Sigma_s(l,g_s_index-1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_s_index-1,l,m,:))

                        end do

                        do g_s_index = g_index,(N%Group-1) ! Downscatter

                            Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + REAL((2.0_8*l + 1),8)*Y_l_m*(Properties%Elements(i)%Sigma_s(l,g_s_index+1,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Legendre_Flux(g_s_index+1,l,m,:))

                        end do

                    end do

                end do

            end do

        end do

        end if

    end subroutine Anisotropic_Scattering

    subroutine Up_Wind_Source(Properties, g_index, ang, i)

        type(PropertiesType), intent(inout) :: Properties

        integer, intent(in) :: i, g_index, ang
        integer             :: side_index, neighbour_element

        Properties%Elements(i)%Sweep_Source(g_index,ang,:) = 0.0_8

        do side_index = 1,Properties%Elements(i)%Number_of_Sides

            neighbour_element = Properties%Elements(i)%Neighbours(side_index,1)

            if (neighbour_element == 0) then

                Properties%Elements(i)%Sweep_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + matmul(Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Properties%Elements(i)%Sides(side_index)%Boundary(g_index,ang,:))

            else

                Properties%Elements(i)%Sweep_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + matmul(Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Properties%Elements(neighbour_element)%Flux(g_index,ang,:))

            end if

        end do

        Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + Properties%Elements(i)%Source(g_index,ang,:)

    end subroutine Up_Wind_Source

    subroutine Calculate_Scalar_Flux(Properties, N)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang, w_index

        real(kind = 8), dimension(:), allocatable :: w

        if (N%D == 1) allocate(w(N%Ordinates))
        if (N%D == 2) allocate(w(N%Ordinates/4))
        if (N%D == 3) allocate(w(N%Ordinates/8))

        if (N%D == 1) call Calculate_w_1D(w)
        if (N%D == 2) call Calculate_w(w)
        if (N%D == 3) call Calculate_w(w)

        do i = 1,N%Element

            Properties%Elements(i)%Scalar_Flux = 0.0_8

            do g_index = 1,N%Group

                do ang = 1,N%Ordinates

                    if (N%D == 1) w_index = ang

                    if (N%D == 2) w_index = MOD(ang,N%Ordinates/4)
                    if (N%D == 3) w_index = MOD(ang,N%Ordinates/8)

                    if (w_index == 0) then

                        if (N%D == 2) w_index = N%Ordinates/4
                        if (N%D == 3) w_index = N%Ordinates/8

                    end if

                    Properties%Elements(i)%Scalar_Flux(g_index,:) = Properties%Elements(i)%Scalar_Flux(g_index,:) + (1.0_8/(2.0_8**N%D))*w(w_index)*Properties%Elements(i)%Flux(g_index,ang,:)

                end do

            end do

        end do

    end subroutine Calculate_Scalar_Flux

    subroutine Input_Results(Results, Properties, N)

        type(ResultsType), intent(inout)    :: Results
        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, dimension(N%N) :: Counter

        integer :: i, g_index, node_index

        Counter = 0

        Results%Scalar_Flux = 0.0_8

        Results%Angular_Flux = 0.0_8

        do i = 1,N%Element

            do node_index = 1,Properties%Elements(i)%Number_of_Nodes

                do g_index = 1,N%Group

                    Results%Scalar_Flux(g_index,Properties%Elements(i)%Cell_Pointers(node_index)) = Results%Scalar_Flux(g_index,Properties%Elements(i)%Cell_Pointers(node_index)) + Properties%Elements(i)%Scalar_Flux(g_index,node_index)

                    Results%Angular_Flux(g_index,Properties%Elements(i)%Cell_Pointers(node_index),:) = Results%Angular_Flux(g_index,Properties%Elements(i)%Cell_Pointers(node_index),:) + Properties%Elements(i)%Flux(g_index,:,node_index)

                end do

                Counter(Properties%Elements(i)%Cell_Pointers(node_index)) = Counter(Properties%Elements(i)%Cell_Pointers(node_index)) + 1

            end do

        end do

        do i = 1,N%N

            Results%Scalar_Flux(:,i) = Results%Scalar_Flux(:,i)/Counter(i)

            Results%Angular_Flux(:,i,:) = Results%Angular_Flux(:,i,:)/Counter(i)

        end do

    end subroutine Input_Results

    subroutine Calculate_Total_Source(Total_Source, N, Properties)

        type(PropertiesType), intent(in) :: Properties
        type(NType), intent(in)          :: N

        real(kind = 8),    intent(inout) :: Total_Source

        integer :: i

        Total_Source = 0.0_8

        do i = 1,N%Element

            Total_Source = Total_Source + SUM(Properties%Elements(i)%Source)

        end do

    end subroutine Calculate_Total_Source

    subroutine Create_Matrices(Properties, N)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang, mu_counter, p, sign

        integer, dimension(N%Ordinates) :: p_indices, q_indices, w_indices, mu_indices

        real(kind = 8), dimension(N%Angle/2) :: mu
        real(kind = 8), dimension(N%Angle)   :: mu_1D

        real(kind = 8), dimension(N%Angle)  :: xsi

        real(kind = 8) :: mu_val, eta_val, xi_val

        do i = 1,N%Element

            if(N%D == 1) call Calculate_Jacobian_1D(Properties, N%Degree, i)
            if(N%D == 2) call Calculate_Jacobian_2D(Properties, N%Degree, i)
            if(N%D == 3) call Calculate_Jacobian_3D(Properties, N%Degree, i)

            if(N%D == 1) call Construct_Streaming_and_Mass_Matrix(Properties, N, i)
            if(N%D == 2) call Construct_Mass_Matrix_2D(Properties, N, i)
            if(N%D == 3) call Construct_Mass_Matrix_3D(Properties, N, i)

            if (Properties%Case == 0) then

                if (Properties%Elements(i)%Cell_Type== 3 .or. Properties%Elements(i)%Cell_Type== 21) then

                    call Calculate_1D_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                else if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    call Calculate_Tri_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                    if (Properties%g == 1) Properties%Elements(i)%Source_Vector = Properties%Elements(i)%Source_Vector*(Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1) + Properties%Elements(i)%Coordinates(3,1))/3.0_8

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                    call Calculate_Quad_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                    if (Properties%g == 1) Properties%Elements(i)%Source_Vector = Properties%Elements(i)%Source_Vector*(Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1) + Properties%Elements(i)%Coordinates(3,1) + Properties%Elements(i)%Coordinates(4,1))/4.0_8

                else if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    call Calculate_Hex_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

                    call Calculate_Tet_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    call Calculate_Pris_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                else if (Properties%Elements(i)%Cell_Type == 14) then

                    call Calculate_Pyr_Source_Vector(Properties,N,i,Properties%Elements(i)%Source_Vector)

                end if

            end if

        end do

        if (N%D == 3) then
        call Calculate_mu(mu)
        !$OMP PARALLEL DO PRIVATE(ang,g_index,i,mu_val,eta_val,xi_val)
        do ang = 1,N%Ordinates

            call Ordinates(ang, N%Angle, mu, mu_val, eta_val, xi_val)

            do g_index = 1,N%Group

                do i = 1,N%Element

                    call Construct_Total_Matrix_3D(Properties, N, i, g_index, ang, mu_val, eta_val, xi_val)

                end do

            end do

        end do
        !$OMP END PARALLEL DO
        else if (N%D == 2) then
        call Calculate_mu(mu)
        if (Properties%g == 0) then
        !$OMP PARALLEL DO PRIVATE(ang,g_index,i,mu_val,eta_val)
        do ang = 1,N%Ordinates

            call Ordinates(ang, N%Angle, mu, mu_val, eta_val, xi_val)

            do g_index = 1,N%Group

                do i = 1,N%Element

                    call Construct_Total_Matrix_2D(Properties, N, i, g_index, ang, mu_val, eta_val)

                end do

            end do

        end do
        !$OMP END PARALLEL DO
        else if (Properties%g == 1) then
        call Create_Ordinate_Arrays_RZ(N, p_indices, q_indices, w_indices, mu_indices)
        call calculate_xsi(xsi)
        !$OMP PARALLEL DO PRIVATE(ang,g_index,i,mu_counter,p,sign,mu_val,eta_val)
        do ang = 1,N%Ordinates

            if (ang == 1) sign = -1
            if (ang == N%Ordinates/4 + 1) sign = 1
            if (ang == N%Ordinates/2 + 1) sign = -1
            if (ang == 3*N%Ordinates/4 + 1) sign = 1

            p = p_indices(ang)
            mu_counter = mu_indices(ang)

            mu_val = sign*mu(mu_counter)
            eta_val = xsi(p)

            do g_index = 1,N%Group

                do i = 1,N%Element

                    call Construct_Total_Matrix_2D(Properties, N, i, g_index, ang, mu_val, eta_val)

                end do

            end do

        end do
        !$OMP END PARALLEL DO
        end if
        else if (N%D == 1) then
        call Calculate_mu_1D(mu_1D)
        !$OMP PARALLEL DO PRIVATE(ang,g_index,i,mu_val)
        do ang = 1,N%Ordinates

            mu_val = mu_1D(ang)

            do g_index = 1,N%Group

                do i = 1,N%Element

                    call Construct_Total_Matrix_1D(Properties, i, g_index, ang, mu_val)

                end do

            end do

        end do
        end if

    end subroutine Create_Matrices

    subroutine Calculate_Boundaries(Properties, N)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang, side_index, neighbour_element

        do g_index = 1, N%Group

            do ang = 1, N%Ordinates

                do i = 1, N%Element

                    do side_index = 1, Properties%Elements(i)%Number_of_Sides

                        neighbour_element = Properties%Elements(i)%Neighbours(side_index,1)

                        if (neighbour_element == 0) then

                            if (N%D == 1) call Boundary_Conditions_1D(Properties, N, i, side_index, g_index, ang, Properties%Elements(i)%Sides(side_index)%Boundary(g_index,ang,:))
                            if (N%D == 2) call Boundary_Conditions_2D(Properties, N, i, side_index, g_index, ang, Properties%Elements(i)%Sides(side_index)%Boundary(g_index,ang,:))
                            if (N%D == 3) call Boundary_Conditions_3D(Properties, N, i, side_index, g_index, ang, Properties%Elements(i)%Sides(side_index)%Boundary(g_index,ang,:))

                        else

                            Properties%Elements(i)%Sides(side_index)%Boundary(g_index,ang,:) = 0.0_8

                        end if

                    end do

                end do

            end do

        end do

    end subroutine Calculate_Boundaries

end module m_Solver