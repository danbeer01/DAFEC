module m_RZ
!
! Purpose:
! To solve the transport equation in 2D cylindrical (RZ) coordinates using the Finite Element Method.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
    use m_Gauss_Points
    use m_Read_Properties
    use m_Calculate_mu_w
    use m_Create_Quadrilateral_Shape_Functions
    use m_Create_Triangular_Shape_Functions
    use m_Small_Matrix_Solver

    implicit none

    contains

    subroutine RZ_Solver(Properties, N, Sweep_Order, lambda_new)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        integer, dimension(:,:), intent(in)   :: Sweep_Order

        real(8), dimension(N%Angle,N%Angle+1) :: alpha
        integer                               :: i, g_index, ang, iter, element_index
        integer                               :: p, q, mu_counter, w_counter

        real(kind = 8), dimension(N%Angle/2)  :: mu
        real(kind = 8), dimension(N%Angle)    :: xsi
        real(kind = 8), dimension(N%Ordinates/4) :: w

        integer, dimension(N%Ordinates) :: p_indices, q_indices, w_indices, mu_indices

        real(kind = 8), dimension(N%Element,4,4) :: L

        real(kind = 8), dimension(N%Group,N%Angle,N%Element,4) :: Half_Flux

        real(kind = 8)                      :: lambda_old = 1.1_8, lambda_new

        real(kind = 8)                      :: Total_Source, Total_Source_new

        integer                             :: max_iter = 10000 ! maximum number of iterations
        real(kind = 8)                      :: tol = 1.0e-6 ! tolerance for convergence

        write(*, '(A)') "*******************************"

        write(*, '(A)') "*                             *"

        write(*, '(A, A, A)') "*  ", "WARNING: RZ CODE IS SHIT ", "  *"

        write(*, '(A)') "*                             *"

        write(*, '(A)') "*******************************"
        
        call Calculate_mu(mu)

        call Calculate_w(w)

        call Calculate_xsi(xsi)

        call Calculate_alpha(N, alpha, w, xsi)

        call Create_Ordinate_Arrays_RZ(N, p_indices, q_indices, w_indices, mu_indices)

        call Construct_RZ_Matrix(Properties, N, alpha, L, mu, w, p_indices, q_indices, w_indices, mu_indices)

        iter = 0

        do while (ABS((lambda_new - lambda_old)/(lambda_old)) > tol)

            lambda_old = lambda_new

            call Calculate_Half_Flux(Properties, N, Half_Flux, Sweep_Order(N%Ordinates+1:N%Ordinates+N%Angle,:), xsi, lambda_new)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_RZ_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_RZ_Total_Source(Total_Source, N, Properties)

            do ang = 1, N%Ordinates

                if (ang == 1) mu = -mu
                if (ang == N%Ordinates/4 + 1) mu = -mu
                if (ang == N%Ordinates/2 + 1) mu = -mu
                if (ang == 3*N%Ordinates/4 + 1) mu = -mu

                w_counter = w_indices(ang)
                p = p_indices(ang)
                q = q_indices(ang)
                mu_counter = mu_indices(ang)

                do g_index = 1, N%Group

                    do i = 1, N%Element

                        element_index = Sweep_Order(ang,i)

                        call Up_Wind_RZ_Source(Properties, N, g_index, ang, element_index, p, q, w(w_counter), alpha, Half_Flux(g_index,p,element_index,:), L(element_index,:,:))

                        call Solve_Matrix(Properties%Elements(element_index)%K_Matrix(g_index,ang,:,:), Properties%Elements(element_index)%Total_Source(g_index,ang,:), Properties%Elements(element_index)%Flux(g_index,ang,:))

                        Half_Flux(g_index,p,element_index,:) = 2.0_8*Properties%Elements(element_index)%Flux(g_index,ang,:) - Half_Flux(g_index,p,element_index,:)

                    end do

                end do

            end do

            call Calculate_RZ_Scalar_Flux(Properties, N, w_indices)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_RZ_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_RZ_Total_Source(Total_Source_new, N, Properties)

            lambda_new=lambda_old*(Total_Source_new/Total_Source) ! Calculate the new eigenvalue

            iter = iter + 1

            if (iter > max_iter) then

                write(*,*) 'Maximum number of iterations reached'
                stop

            end if

        end do

        print *, 'Iterations', iter

    end subroutine RZ_Solver

    subroutine Calculate_RZ_Total_Source(Total_Source, N, Properties)

        type(PropertiesType), intent(in) :: Properties
        type(NType), intent(in)          :: N

        real(kind = 8),    intent(inout) :: Total_Source

        integer :: i

        Total_Source = 0.0_8

        do i = 1,N%Element

            Total_Source = Total_Source + SUM(Properties%Elements(i)%Source)

        end do

    end subroutine Calculate_RZ_Total_Source

    subroutine Calculate_RZ_Scalar_Flux(Properties, N, w_indices)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang

        integer, dimension(:) :: w_indices

        real(kind = 8), dimension(N%Ordinates/4) :: w

        call calculate_w(w)

        do i = 1, N%Element

            Properties%Elements(i)%Scalar_Flux = 0.0_8

            do g_index = 1, N%Group

                do ang = 1, N%Ordinates

                    Properties%Elements(i)%Scalar_Flux(g_index,:) = Properties%Elements(i)%Scalar_Flux(g_index,:) + 0.25_8*w(w_indices(ang))*Properties%Elements(i)%Flux(g_index,ang,:)

                end do

            end do

        end do

    end subroutine Calculate_RZ_Scalar_Flux

    subroutine Calculate_RZ_Source(Properties, N, k_eff, i)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        real(kind = 8) :: k_eff

        integer, intent(in) :: i
        integer :: ang, g_index, g_s_index

        do ang = 1, N%Ordinates

            do g_index = 1, N%Group

                if (Properties%Case == 0) then

                    Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Sigma_f(g_index)*Properties%Elements(i)%Source_Vector(:) + Properties%Elements(i)%Sigma_s(0,g_index,g_index)*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_index,:))

                else if (Properties%Case == 1) then

                    Properties%Elements(i)%Source(g_index,ang,:) = (Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index) + Properties%Elements(i)%Sigma_s(0,g_index,g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_index,:))

                end if

                do g_s_index = g_index,2,-1 ! Upscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*Properties%Elements(i)%Scalar_Flux(g_s_index-1,:)*Properties%Elements(i)%Source_Vector(:)

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*Properties%Elements(i)%Scalar_Flux(g_s_index-1,:)*Properties%Elements(i)%Source_Vector(:)

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index-1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_s_index-1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                    end if

                end do

                do g_s_index = g_index,(N%Group-1) ! Downscatter

                    if (Properties%Case == 0) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*Properties%Elements(i)%Scalar_Flux(g_s_index+1,:)*Properties%Elements(i)%Source_Vector(:)

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*Properties%Elements(i)%Scalar_Flux(g_s_index+1,:)*Properties%Elements(i)%Source_Vector(:)

                    else if (Properties%Case == 1) then

                        if (Properties%Adjoint == 0) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index+1))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                        if (Properties%Adjoint == 1) Properties%Elements(i)%Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_s_index+1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                    end if

                end do

            end do

        end do

    end subroutine Calculate_RZ_Source

    subroutine Calculate_Half_Flux(Properties, N, Half_Flux, Sweep_Order, xsi, k_eff)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        integer, dimension(:,:) :: Sweep_Order

        real(kind = 8), dimension(:,:,:,:), intent(inout) :: Half_Flux

        real(kind = 8), dimension(:), intent(in) :: xsi

        real(kind = 8), dimension(4,4) :: T, K, M, F_out

        real(kind = 8), dimension(4) :: S, Source

        integer :: i, g_index, p, g_s_index, element_index, side_index, neighbour_element

        real(kind = 8) :: k_eff

        real(kind = 8), dimension(:), allocatable :: Boundary

        Half_Flux = 0.0_8

        do p = 1, N%Angle

            do g_index = 1, N%Group
                
                do i = 1,N%Element

                    element_index = Sweep_Order(p,i)

                    call Construct_RZ_F_out_Matrix(Properties, N, element_index, -sqrt(1 - xsi(p)**2), xsi(p), F_out)

                    call Construct_RZ_F_in_Matrix(Properties, N, element_index, N%Ordinates + p, -sqrt(1 - xsi(p)**2), xsi(p))

                    M = Create_M(Properties,N,element_index)

                    K = Create_K(Properties,N,element_index,xsi(p))

                    S = Create_S(Properties,N,element_index)

                    T = - K + Properties%Elements(element_index)%Sigma_t(g_index)*M + F_out

                    if (Properties%Case == 0) then

                        Source = Properties%Elements(element_index)%Sigma_f(g_index)*S + Properties%Elements(element_index)%Sigma_s(0,g_index,g_index)*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_index,:))

                    else if (Properties%Case == 1) then

                        Source = (Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(element_index)%Sigma_f(g_index) + Properties%Elements(element_index)%Sigma_s(0,g_index,g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_index,:))

                    end if

                    do g_s_index = g_index,2,-1 ! Upscatter

                        if (Properties%Case == 0) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index-1,g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index-1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index-1,g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index-1,:))

                        else if (Properties%Case == 1) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(element_index)%Sigma_f(g_s_index-1))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index-1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_s_index-1)*(1.0_8/k_eff)*Properties%Elements(element_index)%Sigma_f(g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index-1,:))

                        end if

                    end do

                    allocate(Boundary(Properties%Elements(i)%Number_of_Nodes))

                    do g_s_index = g_index,(N%Group-1) ! Downscatter

                        if (Properties%Case == 0) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index+1,g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index+1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index+1,g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index+1,:))

                        else if (Properties%Case == 1) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(element_index)%Sigma_f(g_s_index+1))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index+1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(element_index)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_s_index+1)*(1.0_8/k_eff)*Properties%Elements(element_index)%Sigma_f(g_index))*matmul(M,Properties%Elements(element_index)%Scalar_Flux(g_s_index+1,:))

                        end if

                    end do

                    do side_index = 1,Properties%Elements(element_index)%Number_of_Sides

                        neighbour_element = Properties%Elements(element_index)%Neighbours(side_index,1)
            
                        if (neighbour_element == 0) then
            
                            call RZ_Boundary_Conditions(Properties, N, element_index, side_index, g_index, p, 1, Boundary)
            
                            Source = Source + matmul(Properties%Elements(element_index)%Sides(side_index)%F_in_Matrix(N%Ordinates+p,:,:),Boundary)
            
                        else
            
                            Source = Source + matmul(Properties%Elements(element_index)%Sides(side_index)%F_in_Matrix(N%Ordinates+p,:,:),Half_Flux(g_index,p,neighbour_element,:))
            
                        end if
            
                    end do

                    deallocate(Boundary)

                    call Solve_Matrix(T, Source, Half_Flux(g_index,p,element_index,:))

                end do

            end do

        end do

    end subroutine Calculate_Half_Flux

    subroutine Construct_RZ_Matrix(Properties, N, alpha, L, mu, w, p_indices, q_indices, w_indices, mu_indices)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, dimension(:), intent(in) :: p_indices, q_indices, w_indices, mu_indices

        integer :: i, g_index, ang, q, mu_counter, p, w_counter

        real(kind = 8), dimension(:,:), intent(in) :: alpha

        real(kind = 8), dimension(:,:,:) :: L

        real(kind = 8), dimension(:), intent(inout) :: mu, w

        do i = 1, N%Element

            L(i,:,:) = Create_M(Properties,N,i)

            do g_index = 1, N%Group

                do ang = 1, N%Ordinates

                    if (ang == 1) mu = -mu
                    if (ang == N%Ordinates/4 + 1) mu = -mu
                    if (ang == N%Ordinates/2 + 1) mu = -mu
                    if (ang == 3*N%Ordinates/4 + 1) mu = -mu

                    p = p_indices(ang)
                    q = q_indices(ang)
                    mu_counter = mu_indices(ang)
                    w_counter = w_indices(ang)

                    Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = Properties%Elements(i)%K_Matrix(g_index,ang,:,:) + (mu(mu_counter) + (2.0_8/w(w_counter))*alpha(p,q+1))*L(i,:,:)

                end do

            end do

        end do

        mu = abs(mu)

    end subroutine Construct_RZ_Matrix

    subroutine Up_Wind_RZ_Source(Properties, N, g_index, ang, i, p, q, w, alpha, Half_Flux, M)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                 :: i, ang, g_index, p, q

        real(kind = 8), intent(in)          :: w

        integer                             :: side_index, neighbour_element

        real(kind = 8), dimension(:,:), intent(in) :: M

        real(kind = 8), dimension(:,:), intent(in) :: alpha

        real(kind = 8), dimension(:) :: Half_Flux

        real(kind = 8), dimension(:), allocatable :: Boundary

        real(kind = 8) :: r

        Properties%Elements(i)%Sweep_Source(g_index,ang,:) = 0.0_8

        allocate(Boundary(Properties%Elements(i)%Number_of_Nodes))

        Boundary = 0.0_8

        do side_index = 1,Properties%Elements(i)%Number_of_Sides

            ! r = (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(side_index,1),1) + Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(side_index,2),1))/2.0_8
            r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1) + Properties%Elements(i)%Coordinates(3,1) + Properties%Elements(i)%Coordinates(4,1))/4.0_8

            neighbour_element = Properties%Elements(i)%Neighbours(side_index,1)

            if (neighbour_element == 0) then

                call RZ_Boundary_Conditions(Properties, N, i, side_index, g_index, p, q, Boundary)

                Properties%Elements(i)%Sweep_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + matmul(Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Boundary)

            else

                Properties%Elements(i)%Sweep_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + r*matmul(Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:),Properties%Elements(neighbour_element)%Flux(g_index,ang,:))

            end if

        end do

        Properties%Elements(i)%Sweep_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + (alpha(p,q+1) + alpha(p,q))*(1.0_8/w)*matmul(M,Half_Flux)

        deallocate(Boundary)

        Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Sweep_Source(g_index,ang,:) + Properties%Elements(i)%Source(g_index,ang,:)

    end subroutine Up_Wind_RZ_Source

    subroutine RZ_Boundary_Conditions(Properties, N, element_index, side_index, g_index, p, q, Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: element_index, g_index, p, q, side_index

        integer :: p_new, q_new, reflect_ang, N_q_Level_p

        integer :: p_eff, q_eff

        integer :: i, j, k

        real(kind = 8), dimension(:) :: Boundary

        if (Properties%Elements(element_index)%Neighbours(side_index,2) == 2) then

        if (p <= N%Angle/2) then

            N_q_Level_p = 2*p

        else if (p > N%Angle/2) then

            N_q_Level_p = 2*((N%Angle+1)-p)

        end if

        p_new = p

        q_new = N_q_Level_p - (q-1)

        if (p <= N%Angle/2) p_new = N%Angle + 1 - p

        p_eff = p_new - N%Angle/2

        q_eff = q_new - N_q_Level_p/2

        reflect_ang = 0
        i = 1
        j = 1

        do
            reflect_ang = reflect_ang + 1

            if (i == p_eff .and. j == q_eff) exit

            j = j + 1

            k = (N%Angle/2+1-i)

            if (j > k) then

                j = 1
                i = i + 1

            end if

            if(i > N%Angle/2) exit

        end do

        if (p <= N%Angle/2) then

            reflect_ang = reflect_ang + 3*(N%Angle+2)*N%Angle/8

        else if (p > N%Angle/2) then

            reflect_ang = reflect_ang + (N%Angle+2)*N%Angle/8

        end if

        if (Properties%RBC == 1) then

            Boundary = Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%RBC == 2) then

            Boundary = 0.0_8

        else if (Properties%RBC == 3) then

            Boundary = Properties%alpha*Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%RBC == 4) then

            print *, 'Periodic Boundary Conditions not possible'
            stop

        else if (Properties%RBC == 5) then

            Boundary = Properties%Q_s

        end if

        else if (Properties%Elements(element_index)%Neighbours(side_index,2) == 3) then

        p_new = N%Angle + 1 - p

        q_new = q

        p_eff = p_new - N%Angle/2

        N_q_Level_p = 2*p

        if (q <= N%Angle/2) then

            q_eff = q_new

        else if (q > N%Angle/2) then

            q_eff = q_new - N_q_Level_p/2

        end if

        reflect_ang = 0
        i = 1
        j = 1

        do
            reflect_ang = reflect_ang + 1

            if (i == p_eff .and. j == q_eff) exit

            j = j + 1

            k = (N%Angle/2+1-i)

            if (j > k) then

                j = 1
                i = i + 1

            end if

            if(i > N%Angle/2) exit

        end do

        if (q > N%Angle/2) then

            reflect_ang = reflect_ang + (N%Angle+2)*N%Angle/8

        end if

        if (Properties%TBC == 1) then

            Boundary = Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%TBC == 2) then

            Boundary = 0.0_8

        else if (Properties%TBC == 3) then

            Boundary = Properties%alpha*Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%TBC == 4) then

            print *, 'Periodic Boundary Conditions not possible'
            stop

        else if (Properties%TBC == 5) then

            Boundary = Properties%Q_s

        end if

        else if (Properties%Elements(element_index)%Neighbours(side_index,2) == 4) then

        p_new = N%Angle + 1 - p

        q_new = q

        p_eff = N%Angle/2 + 1 - p_new

        N_q_Level_p = 2*((N%Angle+1)-p)

        if (q <= N%Angle/2) then

            q_eff = q_new

        else if (q > N%Angle/2) then

            q_eff = q_new - N_q_Level_p/2

        end if

        reflect_ang = 0
        i = 1
        j = 1

        do
            reflect_ang = reflect_ang + 1

            if (i == p_eff .and. j == q_eff) exit

            j = j + 1

            k = (N%Angle/2+1-i)

            if (j > k) then

                j = 1
                i = i + 1

            end if

            if(i > N%Angle/2) exit

        end do

        if (q > N%Angle/2) then

            reflect_ang = reflect_ang + 3*(N%Angle+2)*N%Angle/8

        else if (q <= N%Angle/2) then

            reflect_ang = reflect_ang + 2*(N%Angle+2)*N%Angle/8

        end if

        if (Properties%BBC == 1) then

            Boundary = Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%BBC == 2) then

            Boundary = 0.0_8

        else if (Properties%BBC == 3) then

            Boundary = Properties%alpha*Properties%Elements(element_index)%Flux(g_index,reflect_ang,:)

        else if (Properties%BBC == 4) then

            print *, 'Periodic Boundary Conditions not possible'
            stop

        else if (Properties%BBC == 5) then

            Boundary = Properties%Q_s

        end if

        end if

    end subroutine RZ_Boundary_Conditions

    subroutine Calculate_alpha(N, alpha, w, mu)

        type(NType), intent(in)  :: N

        real(kind = 8), dimension(:,:), intent(inout) :: alpha

        real(kind = 8), dimension(:), intent(in) :: w, mu

        integer :: p, q, Num_q_Level_p, w_counter_1, w_counter_2, w_counter

        alpha = 0.0_8

        Num_q_Level_p = 1
  
        do p = 1, N%Angle
  
          do q = 2, 2*Num_q_Level_p + 1
  
            if (q - 1 <= Num_q_Level_p) then

                if (p <= N%Angle/2) then

                    if (p == 1) then
      
                        w_counter_1 = N%Ordinates/4
          
                    else 
          
                        w_counter_1 = w_counter_1 - 1
          
                    end if
      
                else
  
                    if (p == N%Angle/2 + 1 .and. q-1 == 1) then
        
                        w_counter_1 = Num_q_Level_p
                        w_counter = w_counter_1
        
                    else if (q-1 == 1) then
            
                        w_counter_1 = w_counter + Num_q_Level_p
                        w_counter = w_counter_1
        
                    else
            
                        w_counter_1 = w_counter_1 - 1
            
                    end if
        
                end if

                alpha(p,q) = alpha(p,q-1) - w(w_counter_1)*(mu(q+N%Angle/2-Num_q_Level_p-1))
  
            else

                if (p <= N%Angle/2) then

                    if (p == 1) then
      
                        w_counter_2 = N%Ordinates/4
                        w_counter = N%Ordinates/4
              
                    else if (q-1 == Num_q_Level_p + 1) then
            
                        w_counter_2 = w_counter - Num_q_Level_p
                        w_counter = w_counter_2
          
                    else 
          
                        w_counter_2 = w_counter_2 + 1
          
                    end if
      
                else
        
                    if (p == N%Angle/2 + 1 .and. q-1 == N%Angle/2 + 1) then
        
                        w_counter_2 = 1
            
                    else 
            
                        w_counter_2 = w_counter_2 + 1
            
                    end if

                end if
        
              alpha(p,q) = alpha(p,q-1) - w(w_counter_2)*(mu(q+N%Angle/2-Num_q_Level_p-1))
  
            end if
  
          end do

          if (p < N%Angle/2) then

            Num_q_Level_p = Num_q_Level_p + 1
  
          else if (p > N%Angle/2) then
  
            Num_q_Level_p = Num_q_Level_p - 1
  
          end if
  
        end do

    end subroutine Calculate_alpha

    Function Create_K(Properties,N,i,xsi) Result(K)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in)  :: i

        real(kind = 8), intent(in) :: xsi

        integer :: j, index_1, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), allocatable :: K

        real(kind = 8), dimension(:,:,:), allocatable :: dSFMat, dSFMatT

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then
            if(N%Degree == 1) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if
        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
            Num_Gauss_Points = Num_Nodes
        end if

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(K(Num_Nodes,Num_Nodes))

        allocate(dSFMatT(Num_Gauss_Points,Num_Nodes,2), dSFMat(Num_Gauss_Points,2,Num_Nodes))

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        K = 0.0_8

        do j = 1, Num_Gauss_Points

            if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                do index_1 = 1, Num_Nodes

                    dSFMat(j,1,index_1) = -SQRT(1 - xsi**2)*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))
                    dSFMat(j,2,index_1) = xsi*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))
    
                    dSFMatT(j,index_1,1) = Generate_Tri_Shape_Functions_Derivative_xi(N%Degree, index_1, eta(j), xi(j))
                    dSFMatT(j,index_1,2) = Generate_Tri_Shape_Functions_Derivative_eta(N%Degree, index_1, eta(j), xi(j))

                end do

                K = K + 0.5_8*w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

            else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                do index_1 = 1, Num_Nodes

                    dSFMat(j,1,index_1) = -SQRT(1 - xsi**2)*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j))
                    dSFMat(j,2,index_1) = xsi*Generate_Quadrilateral_Shape_Functions(Properties, N%Degree, index_1, eta(j), xi(j))

                    dSFMatT(j,index_1,1) = Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, index_1, eta(j), xi(j))
                    dSFMatT(j,index_1,2) = Generate_Shape_Functions_Derivative_eta(Properties, N%Degree, index_1, eta(j), xi(j))

                end do

                K = K + w(j)*ABS(Properties%Elements(i)%Det_Jacobian(j))*matmul(matmul(dSFMatT(j,:,:), transpose(Properties%Elements(i)%Inverse_Jacobian(j,:,:))), dSFMat(j,:,:))

            end if

        end do

    end Function Create_K

    Function Create_M(Properties,N,i) Result(M)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i

        integer             :: j, index_1, index_2, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:,:), allocatable :: M

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then
            if(N%Degree == 1) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if
        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
            Num_Gauss_Points = Num_Nodes
        end if

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(M(Num_Nodes,Num_Nodes))

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        M = 0.0_8

        do j = 1, Num_Gauss_Points

            do index_1 = 1, Num_Nodes

                do index_2 = 1, Num_Nodes

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) M(index_1,index_2) = M(index_1,index_2) + 0.5_8*w(j)*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))*Generate_Triangular_Shape_Functions(N%Degree, index_2, eta(j), xi(j))*Properties%Elements(i)%Det_Jacobian(j)
                    if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) M(index_1,index_2) = M(index_1,index_2) + w(j)*Generate_Quadrilateral_Shape_Functions(Properties,N%Degree, index_1, eta(j), xi(j))*Generate_Quadrilateral_Shape_Functions(Properties,N%Degree, index_2, eta(j), xi(j))*Properties%Elements(i)%Det_Jacobian(j)

                end do

            end do

        end do

    end Function Create_M

    Function Create_S(Properties,N,i) Result(S)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in)  :: i

        integer             :: j, index_1, Num_Gauss_Points, Num_Nodes

        real(kind = 8), dimension(:), allocatable :: S

        real(kind = 8), dimension(:), allocatable :: xi, eta, w

        Num_Nodes = Properties%Elements(i)%Number_of_Nodes

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then
            if(N%Degree == 1) then
                Num_Gauss_Points = 3
            else
                Num_Gauss_Points = 7
            end if
        else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
            Num_Gauss_Points = Num_Nodes
        end if

        allocate(xi(Num_Gauss_Points), eta(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(S(Num_Nodes))

        if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) call Generate_Triangular_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) call Generate_2D_Quad_Gauss_Points(Num_Gauss_Points, eta, xi, w)

        S = 0.0_8

        do j = 1, Num_Gauss_Points

            do index_1 = 1, Num_Nodes

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) S(index_1) = S(index_1) + 0.5_8*w(j)*Generate_Triangular_Shape_Functions(N%Degree, index_1, eta(j), xi(j))*Properties%Elements(i)%Det_Jacobian(j)
                if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) S(index_1) = S(index_1) + w(j)*Generate_Quadrilateral_Shape_Functions(Properties,N%Degree, index_1, eta(j), xi(j))*Properties%Elements(i)%Det_Jacobian(j)

            end do

        end do

    end Function Create_S

    subroutine Create_Ordinate_Arrays_RZ(N, p_indices, q_indices, w_indices, mu_indices)

        type(NType), intent(in)              :: N

        integer, dimension(:) :: p_indices, q_indices, w_indices, mu_indices

        integer :: ang, p, q
        integer :: mu_counter, w_counter, Num_q_Level_p, w_counter_2

        do ang = 1, N%Ordinates

            q = q + 1

            if (ang < N%Ordinates/4 + 1 .or. (ang >= N%Ordinates/2 + 1 .and. ang < 3*N%Ordinates/4 + 1)) then

                mu_counter = mu_counter - 1

            else

                mu_counter = mu_counter + 1

            end if
            
            if (ang == 1) then

                p = N%Angle/2 + 1

                q = 1

                Num_q_Level_p = N%Angle

                w_counter = N%Angle/2

                mu_counter = N%Angle/2

            else if (q == Num_q_Level_p/2 + 1 .or. q == Num_q_Level_p + 1) then

                if (ang < N%Ordinates/2 + 1) then

                    p = p + 1

                else

                    p = p - 1

                end if

                Num_q_Level_p = Num_q_Level_p - 2

                if (ang < N%Ordinates/4 + 1 .or. (ang >= N%Ordinates/2 + 1 .and. ang < 3*N%Ordinates/4 + 1)) then

                    q = 1
                    mu_counter = Num_q_Level_p/2

                else

                    q = Num_q_Level_p/2 + 1
                    mu_counter = 1

                end if

            end if

            if (ang == N%Ordinates/4 + 1) then

                p = N%Angle/2 + 1

                q = N%Angle/2 + 1

                Num_q_Level_p = N%Angle

            else if (ang == N%Ordinates/2 + 1) then

                p = N%Angle/2

                q = 1

                Num_q_Level_p = N%Angle

                mu_counter = N%Angle/2

            else if (ang == 3*N%Ordinates/4 + 1) then

                p = N%Angle/2

                q = N%Angle/2 + 1

                Num_q_Level_p = N%Angle

            end if

            if (q <= Num_q_Level_p/2) then

                if (ang == 1 .or. ang == N%Ordinates/2 + 1) then

                    w_counter = N%Angle/2
                    w_counter_2 = N%Angle/2

                else if (q == 1) then

                    w_counter = w_counter_2 + Num_q_Level_p/2
                    w_counter_2 = w_counter

                else

                    w_counter = w_counter - 1

                end if

                w_indices(ang) = w_counter

            else

                if (ang == N%Ordinates/4 + 1 .or. ang == 3*N%Ordinates/4 + 1) then

                    w_counter = 1

                else 

                    w_counter = w_counter + 1

                end if

                w_indices(ang) = w_counter

            end if

            p_indices(ang) = p

            q_indices(ang) = q

            mu_indices(ang) = mu_counter

        end do

    end subroutine Create_Ordinate_Arrays_RZ

    subroutine Construct_RZ_F_out_Matrix(Properties, N, i, mu, eta, F_out)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i

        real(kind = 8), dimension(:,:) :: F_out

        real(kind = 8) :: Omega_n

        integer :: side_index

        F_out = 0.0_8

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta)

            if (Omega_n > 0) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Tri_Side(Properties,N,i,side_index)

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                    F_out(:,:) = F_out + Omega_n*Integrate_Quad_Side(Properties,N,i,side_index)

                end if

            end if

        end do

    end subroutine Construct_RZ_F_out_Matrix

    subroutine Construct_RZ_F_in_Matrix(Properties, N, i, ang, mu, eta)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind=8), intent(in) :: mu, eta
        integer, intent(in)      :: i, ang
        
        real(kind = 8) :: Omega_n

        integer :: side_index

        do side_index = 1, Properties%Elements(i)%Number_of_Sides

            Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = 0.0_8

            Omega_n = (Properties%Elements(i)%Unit_Vectors(side_index,1)*mu + Properties%Elements(i)%Unit_Vectors(side_index,2)*eta)

            if (Omega_n < 0) then

                if (Properties%Elements(i)%Neighbours(side_index,1) == 0) then

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tri_Side(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Quad_Side(Properties,N,i,side_index)

                    end if

                else

                    if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Tri_Side_F_in(Properties,N,i,side_index)

                    else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then

                        Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:) = abs(Omega_n)*Integrate_Quad_Side_F_in(Properties,N,i,side_index)

                    end if

                end if

            end if

        end do

    end subroutine Construct_RZ_F_in_Matrix

end module m_RZ