module m_Cylindrical
!
! Purpose:
! To solve the transport equation in cylindrical coordinates using the Finite Element Method.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
    use m_Read_Properties
    use m_Calculate_mu_w
    use m_Create_Shape_Functions
    use m_Small_Matrix_Solver

    implicit none

    contains

    subroutine Cylindrical_Solver(Properties, N, lambda_new)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(8), dimension(N%Angle,N%Angle+1) :: alpha
        integer                               :: i, g_index, ang, iter
        integer                               :: p, q, mu_counter, w_counter

        real(kind = 8), dimension(N%Angle/2)  :: mu
        real(kind = 8), dimension(N%Angle)    :: xsi
        real(kind = 8), dimension(N%Ordinates/2) :: w

        integer, dimension(N%Ordinates) :: p_indices, q_indices, w_indices, mu_indices

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Group,N%Angle,N%Element,N%Degree+1) :: Half_Flux

        real(kind = 8)                      :: lambda_old = 1.1_8, lambda_new

        real(kind = 8)                      :: Total_Source, Total_Source_new

        integer                             :: max_iter = 10000 ! maximum number of iterations
        real(kind = 8)                      :: tol = 1.0e-6 ! tolerance for convergence

        call Calculate_mu(mu)

        call Calculate_w(w)

        call Calculate_xsi(xsi)

        call Calculate_alpha(N, alpha, w, xsi)

        call Create_Ordinate_Arrays(N, p_indices, q_indices, w_indices, mu_indices)

        call Construct_Cylindrical_Matrix(Properties, N, alpha, L, mu, w, p_indices, q_indices, w_indices, mu_indices)

        iter = 0

        do while (ABS((lambda_new - lambda_old)/(lambda_old)) > tol)

            lambda_old = lambda_new

            call Calculate_Half_Flux(Properties, N, Half_Flux, xsi, lambda_new)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Cylindrical_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_Cylindrical_Total_Source(Total_Source, N, Properties)

            do g_index = 1, N%Group

                do ang = 1,N%Ordinates

                    if (ang == 1) mu = -mu
                    if (ang == (N%Angle+2)*N%Angle/8 + 1) mu = -mu

                    w_counter = w_indices(ang)
                    p = p_indices(ang)
                    q = q_indices(ang)
                    mu_counter = mu_indices(ang)

                    if (mu(mu_counter) < 0.0_8) then

                        do i = N%Element,1,-1

                            call Up_Wind_Cylindrical_Source(Properties, N, g_index, ang, p, q, i, mu(mu_counter), w(w_counter), alpha, L, Half_Flux(g_index,p,i,:))

                            call Solve_Matrix(Properties%Elements(i)%K_Matrix(g_index,ang,:,:), Properties%Elements(i)%Total_Source(g_index,ang,:), Properties%Elements(i)%Flux(g_index,ang,:))

                            Half_Flux(g_index,p,i,:) = 2.0_8*Properties%Elements(i)%Flux(g_index,ang,:) - Half_Flux(g_index,p,i,:)

                        end do

                    else

                        do i = 1, N%Element

                            call Up_Wind_Cylindrical_Source(Properties, N, g_index, ang, p, q, i, mu(mu_counter), w(w_counter), alpha, L, Half_Flux(g_index,p,i,:))

                            call Solve_Matrix(Properties%Elements(i)%K_Matrix(g_index,ang,:,:), Properties%Elements(i)%Total_Source(g_index,ang,:), Properties%Elements(i)%Flux(g_index,ang,:))

                            Half_Flux(g_index,p,i,:) = 2.0_8*Properties%Elements(i)%Flux(g_index,ang,:) - Half_Flux(g_index,p,i,:)

                        end do

                    end if

                end do

            end do

            call Calculate_Cylindrical_Scalar_Flux(Properties, N, w_indices)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Cylindrical_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_Cylindrical_Total_Source(Total_Source_new, N, Properties)

            lambda_new=lambda_old*(Total_Source_new/Total_Source) ! Calculate the new eigenvalue

            iter = iter + 1

            if (iter > max_iter) then

                write(*,*) 'Maximum number of iterations reached'
                stop

            end if

        end do

        print *, 'Iterations', iter

    end subroutine Cylindrical_Solver

    subroutine Calculate_Cylindrical_Total_Source(Total_Source, N, Properties)

        type(PropertiesType), intent(in) :: Properties
        type(NType), intent(in)          :: N

        real(kind = 8),    intent(inout) :: Total_Source

        integer :: i

        Total_Source = 0.0_8

        do i = 1,N%Element

            Total_Source = Total_Source + SUM(Properties%Elements(i)%Source)

        end do

    end subroutine Calculate_Cylindrical_Total_Source

    subroutine Calculate_Cylindrical_Scalar_Flux(Properties, N, w_indices)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang

        integer, dimension(:) :: w_indices

        real(kind = 8), dimension((N%Angle+2)*N%Angle/8) :: w

        call calculate_w(w)

        do i = 1,N%Element

            Properties%Elements(i)%Scalar_Flux = 0.0_8

            do g_index = 1,N%Group

                do ang = 1,N%Ordinates

                    Properties%Elements(i)%Scalar_Flux(g_index,:) = Properties%Elements(i)%Scalar_Flux(g_index,:) + 0.5_8*w(w_indices(ang))*Properties%Elements(i)%Flux(g_index,ang,:)

                end do

            end do

        end do

    end subroutine Calculate_Cylindrical_Scalar_Flux

    subroutine Calculate_Cylindrical_Source(Properties, N, k_eff, i)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        real(kind = 8) :: k_eff

        integer, intent(in) :: i
        integer :: ang, g_index, g_s_index

        do ang = 1,N%Ordinates

            do g_index = 1,N%Group

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

    end subroutine Calculate_Cylindrical_Source

    subroutine Calculate_Half_Flux(Properties, N, Half_Flux, xsi, k_eff)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind = 8), dimension(:,:,:,:), intent(inout) :: Half_Flux

        real(kind = 8), dimension(:), intent(in) :: xsi

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: T, K, M, F_out, F_in

        real(kind = 8), dimension(N%Degree+1) :: S, Source, Boundary_Flux

        integer :: i, g_index, p, g_s_index

        real(kind = 8) :: dr, k_eff, Boundary

        F_out = 0.0_8

        F_out(1,1) = 1.0_8

        F_in = 0.0_8

        F_in(2,1) = 1.0_8

        Half_Flux = 0.0_8

        Boundary_Flux = 0.0_8

        M = Create_M(Properties,N)

        K = Create_K(Properties,N)

        S = Create_S(Properties,N)

        do p = N%Angle/2+1, N%Angle

            do g_index = 1, N%Group

                call Cylindrical_Boundary_Conditions(Properties, N, g_index, p, 1, Boundary)

                Boundary_Flux(1) = Boundary

                do i = N%Element,1,-1

                    dr = abs(Properties%Elements(i)%Coordinates(1,1) - Properties%Elements(i)%Coordinates(2,1))

                    T = (-2.0_8/dr)*(-F_out - K)*sqrt(1 - xsi(p)**2) + Properties%Elements(i)%Sigma_t(g_index)*M

                    if (Properties%Case == 0) then

                        Source = Properties%Elements(i)%Sigma_f(g_index)*S + Properties%Elements(i)%Sigma_s(0,g_index,g_index)*matmul(M,Properties%Elements(i)%Scalar_Flux(g_index,:))

                    else if (Properties%Case == 1) then

                        Source = (Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index) + Properties%Elements(i)%Sigma_s(0,g_index,g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_index,:))

                    end if

                    do g_s_index = g_index,2,-1 ! Upscatter

                        if (Properties%Case == 0) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                        else if (Properties%Case == 1) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index-1))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index-1,g_index) + Properties%Chi(g_s_index-1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index-1,:))

                        end if

                    end do

                    do g_s_index = g_index,(N%Group-1) ! Downscatter

                        if (Properties%Case == 0) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                        else if (Properties%Case == 1) then

                            if (Properties%Adjoint == 0) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_index)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_s_index+1))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                            if (Properties%Adjoint == 1) Source = Source + (Properties%Elements(i)%Sigma_s(0,g_s_index+1,g_index) + Properties%Chi(g_s_index+1)*(1.0_8/k_eff)*Properties%Elements(i)%Sigma_f(g_index))*matmul(M,Properties%Elements(i)%Scalar_Flux(g_s_index+1,:))

                        end if

                    end do

                    if (i == N%Element) then

                        Source = Source + (2.0_8/dr)*MATMUL(F_in,Boundary_Flux)*sqrt(1 - xsi(p)**2)

                    else

                        Source = Source + (2.0_8/dr)*MATMUL(F_in,Half_Flux(g_index,p,i+1,:))*sqrt(1 - xsi(p)**2)

                    end if

                    call Solve_Matrix(T, Source, Half_Flux(g_index,p,i,:))

                end do

            end do

        end do

    end subroutine Calculate_Half_Flux

    subroutine Construct_Cylindrical_Matrix(Properties, N, alpha, L, mu, w, p_indices, q_indices, w_indices, mu_indices)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, dimension(:), intent(in) :: p_indices, q_indices, w_indices, mu_indices

        integer :: i, g_index, ang, q, mu_counter, p, w_counter

        real(kind = 8), dimension(:,:), intent(in) :: alpha

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: F_out

        real(kind = 8), dimension(:), intent(inout) :: mu, w

        real(kind = 8) :: r_minus, r_plus, dr

        L = Create_L(Properties,N)

        do i = 1, N%Element

            r_minus = Properties%Elements(i)%Coordinates(1,1)

            r_plus = Properties%Elements(i)%Coordinates(2,1)

            dr = abs(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))

            do g_index = 1, N%Group

                do ang = 1, N%Ordinates

                    if (ang == 1) mu = -mu
                    if (ang == (N%Angle+2)*N%Angle/8 + 1) mu = -mu

                    p = p_indices(ang)
                    q = q_indices(ang)
                    mu_counter = mu_indices(ang)
                    w_counter = w_indices(ang)

                    if (mu(mu_counter) < 0.0_8) then

                        F_out = 0.0_8

                        F_out(1,1) = 1.0_8

                        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = (2.0_8*mu(mu_counter)/dr)*(-F_out*r_minus - Properties%Elements(i)%S_Matrix(1,:,:)) + (2.0_8/w(w_counter))*((alpha(p,q+1))*L) + Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix(:,:)

                    else

                        F_out = 0.0_8

                        F_out(2,2) = 1.0_8

                        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = (2.0_8*mu(mu_counter)/dr)*(F_out*r_plus - Properties%Elements(i)%S_Matrix(1,:,:)) + (2.0_8/w(w_counter))*((alpha(p,q+1))*L) + Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix(:,:)

                    end if

                end do

            end do

        end do

        mu = abs(mu)

    end subroutine Construct_Cylindrical_Matrix

    subroutine Up_Wind_Cylindrical_Source(Properties, N, g_index, ang, p, q, i, mu, w, alpha, L, Half_Flux)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                 :: i, ang, g_index, p, q

        real(kind = 8), dimension(:,:), intent(in) :: alpha

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: F_in

        real(kind = 8), intent(in) :: mu, w
        real(kind = 8), dimension(:,:), intent(in) :: L

        real(kind = 8), dimension(N%Degree+1) :: Half_Flux, Boundary_Flux

        real(kind = 8) :: r_minus, r_plus, dr, Boundary

        dr = abs(Properties%Elements(i)%Coordinates(1,1) - Properties%Elements(i)%Coordinates(2,1))

        F_in = 0.0_8

        if (mu < 0.0_8) then

            call Cylindrical_Boundary_Conditions(Properties, N, g_index, p, q, Boundary)

            Boundary_Flux = 0.0_8
    
            Boundary_Flux(1) = Boundary

            F_in(2,1) = 1.0_8

            r_plus = Properties%Elements(i)%Coordinates(2,1)

            if (i == N%Element) then

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(p,q+1) + alpha(p,q))*(1.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) - (2.0_8*mu/dr)*matmul(F_in,Boundary_Flux)*r_plus

            else

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(p,q+1) + alpha(p,q))*(1.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) - (2.0_8*mu/dr)*matmul(F_in,Properties%Elements(i+1)%Flux(g_index,ang,:))*r_plus

            end if

        else

            F_in(1,2) = 1.0_8

            r_minus = Properties%Elements(i)%Coordinates(1,1)

            if (i == 1) then

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(p,q+1) + alpha(p,q))*(1.0_8/w)*MATMUL(L(:,:),Half_Flux(:))

            else

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(p,q+1) + alpha(p,q))*(1.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) + (2.0_8*mu/dr)*matmul(F_in,Properties%Elements(i-1)%Flux(g_index,ang,:))*r_minus

            end if

        end if

    end subroutine Up_Wind_Cylindrical_Source

    subroutine Cylindrical_Boundary_Conditions(Properties, N, g_index, p, q, Right_Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: g_index, p, q

        integer :: p_new, q_new, reflect_ang, N_q_Level_p

        integer :: p_eff, q_eff

        integer :: i, j, k

        real(kind = 8) :: Right_Boundary

        N_q_Level_p = 2*((N%Angle+1)-p)

        p_new = p

        q_new = N_q_Level_p - (q-1)

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

        reflect_ang = reflect_ang + (N%Angle+2)*N%Angle/8

        if (Properties%RBC == 1) then

            Right_Boundary = Properties%Elements(N%Element)%Flux(g_index,reflect_ang,2)

        else if (Properties%RBC == 2) then

            Right_Boundary = 0.0_8

        else if (Properties%RBC == 3) then

            Right_Boundary = Properties%alpha*Properties%Elements(N%Element)%Flux(g_index,reflect_ang,2)

        else if (Properties%RBC == 4) then

            print *, 'Periodic Boundary Conditions not possible'
            stop

        else if (Properties%RBC == 5) then

            Right_Boundary = Properties%Q_s

        end if

    end subroutine Cylindrical_Boundary_Conditions

    subroutine Calculate_alpha(N, alpha, w, mu)

        type(NType), intent(in)  :: N

        real(kind = 8), dimension(:,:), intent(inout) :: alpha

        real(kind = 8), dimension(:), intent(in) :: w, mu

        integer :: p, q, Num_q_Level_p, w_counter_1, w_counter_2, w_counter

        alpha = 0.0_8

        Num_q_Level_p = N%Angle/2
  
        do p = N%Angle/2+1, N%Angle
  
          do q = 2, 2*Num_q_Level_p + 1
  
            if (q - 1 <= Num_q_Level_p) then
  
              if (p == N%Angle/2 + 1 .and. q-1 == 1) then
  
                w_counter_1 = Num_q_Level_p
                w_counter = w_counter_1
  
              else if (q-1 == 1) then
      
                w_counter_1 = w_counter + Num_q_Level_p
                w_counter = w_counter_1
  
              else
    
                w_counter_1 = w_counter_1 - 1
    
              end if
  
              alpha(p,q) = alpha(p,q-1) - w(w_counter_1)*(mu(q+N%Angle/2-Num_q_Level_p-1))
  
            else
  
              if (p == N%Angle/2 + 1 .and. q-1 == N%Angle/2 + 1) then
  
                w_counter_2 = 1
    
              else 
    
                w_counter_2 = w_counter_2 + 1
    
              end if
  
              alpha(p,q) = alpha(p,q-1) - w(w_counter_2)*(mu(q+N%Angle/2-Num_q_Level_p-1))
  
            end if
  
          end do
  
          Num_q_Level_p = Num_q_Level_p - 1
  
        end do

    end subroutine Calculate_alpha

    Function Create_K(Properties,N) Result(K)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: K

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8)                        :: x_val

        K = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        do j = 1, N%Degree+1

            x_val = x(j)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    K(index_1,index_2) = K(index_1,index_2) + w(j)*Generate_Shape_Functions_Derivative(Properties%Isoparametric_Coordinates(:,1), N%Degree, index_1, x_val)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x_val)

                end do

            end do

        end do

    end Function Create_K

    Function Create_M(Properties,N) Result(M)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: M

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8)                        :: x_val

        M = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        do j = 1, N%Degree+1

            x_val = x(j)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    M(index_1,index_2) = M(index_1,index_2) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x_val)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x_val)

                end do

            end do

        end do

    end Function Create_M

    Function Create_S(Properties,N) Result(S)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer             :: j, index_1

        real(kind = 8), dimension(N%Degree+1) :: S

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8)                        :: x_val

        S = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        do j = 1, N%Degree+1

            x_val = x(j)

            do index_1 = 1, N%Degree+1

                S(index_1) = S(index_1) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x_val)

            end do

        end do

    end Function Create_S

    Function Create_L(Properties,N) Result(L)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Degree+1) :: x
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8)                        :: x_val

        L = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        do j = 1, N%Degree+1

            x_val = x(j)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    L(index_1,index_2) = L(index_1,index_2) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x_val)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x_val)

                end do

            end do

        end do

    end Function Create_L

    subroutine Create_Ordinate_Arrays(N, p_indices, q_indices, w_indices, mu_indices)

        type(NType), intent(in)              :: N

        integer, dimension(:) :: p_indices, q_indices, w_indices, mu_indices

        integer :: ang, p, q
        integer :: mu_counter, w_counter, Num_q_Level_p, w_counter_2

        do ang = 1, N%Ordinates

            q = q + 1

            if (ang < (N%Angle+2)*N%Angle/8 + 1) then

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

                p = p + 1

                Num_q_Level_p = Num_q_Level_p - 2

                if (ang < (N%Angle+2)*N%Angle/8 + 1) then

                    q = 1
                    mu_counter = Num_q_Level_p/2

                else

                    q = Num_q_Level_p/2 + 1
                    mu_counter = 1

                end if

            end if

            if (ang == (N%Angle+2)*N%Angle/8 + 1) then

                p = N%Angle/2 + 1

                q = N%Angle/2 + 1

                Num_q_Level_p = N%Angle

            end if

            if (q <= Num_q_Level_p/2) then

                if (ang == 1) then

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

                if (ang == (N%Angle+2)*N%Angle/8 + 1) then

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

    end subroutine Create_Ordinate_Arrays

end module m_Cylindrical