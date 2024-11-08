module m_Spherical
!
! Purpose:
! To solve the transport equation in spherical coordinates using the Finite Element Method.
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

    subroutine Spherical_Solver(Properties, N, lambda_new)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(8), dimension(N%Angle+1) :: alpha
        integer                       :: i, g_index, ang, iter

        real(kind = 8), dimension(N%Angle)  :: mu
        real(kind = 8), dimension(N%Angle)  :: w

        real(kind = 8), dimension(N%Element,N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Group,N%Element,N%Degree+1) :: Half_Flux

        real(kind = 8)                      :: lambda_old = 1.1_8, lambda_new

        real(kind = 8)                      :: Total_Source, Total_Source_new

        integer                             :: max_iter = 10000 ! maximum number of iterations
        real(kind = 8)                      :: tol = 1.0e-8 ! tolerance for convergence

        call Calculate_mu(mu)

        call Calculate_w(w)

        alpha(1) = 0.0_8

        do i = 2, N%Angle + 1

            alpha(i) = alpha(i-1) + mu(i-1)*w(i-1)

        end do

        call Construct_Spherical_Matrix(Properties, N, alpha, L, mu, w)

        iter = 0

        do while (ABS((lambda_new - lambda_old)/(lambda_old)) > tol)

            lambda_old = lambda_new

            call Calculate_Half_Flux(Properties, N, Half_Flux, lambda_new)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Spherical_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_Spherical_Total_Source(Total_Source, N, Properties)

            do g_index = 1, N%Group

                do ang = N%Angle,1,-1

                    if (ang > N%Angle/2) then

                        do i = N%Element,1,-1

                            call Up_Wind_Spherical_Source(Properties, N, g_index, ang, i, mu(ang), w(ang), alpha, L(i,:,:), Half_Flux(g_index,i,:))

                            call Solve_Matrix(Properties%Elements(i)%K_Matrix(g_index,ang,:,:), Properties%Elements(i)%Total_Source(g_index,ang,:), Properties%Elements(i)%Flux(g_index,ang,:))

                            Half_Flux(g_index,i,:) = 2.0_8*Properties%Elements(i)%Flux(g_index,ang,:) - Half_Flux(g_index,i,:)

                        end do

                    else

                        do i = 1, N%Element

                            call Up_Wind_Spherical_Source(Properties, N, g_index, ang, i, mu(ang), w(ang), alpha, L(i,:,:), Half_Flux(g_index,i,:))

                            call Solve_Matrix(Properties%Elements(i)%K_Matrix(g_index,ang,:,:), Properties%Elements(i)%Total_Source(g_index,ang,:), Properties%Elements(i)%Flux(g_index,ang,:))

                            Half_Flux(g_index,i,:) = 2.0_8*Properties%Elements(i)%Flux(g_index,ang,:) - Half_Flux(g_index,i,:)

                        end do

                    end if

                end do

            end do

            call Calculate_Spherical_Scalar_Flux(Properties, N)

            do i = 1,N%Element

                ! Calculate the source terms for each element
                call Calculate_Spherical_Source(Properties, N, lambda_new, i)

            end do

            call Calculate_Spherical_Total_Source(Total_Source_new, N, Properties)

            lambda_new=lambda_old*(Total_Source_new/Total_Source) ! Calculate the new eigenvalue

            iter = iter + 1

            if (iter > max_iter) then

                write(*,*) 'Maximum number of iterations reached'
                stop

            end if

        end do

        print *, 'Iterations', iter

    end subroutine Spherical_Solver

    subroutine Calculate_Spherical_Total_Source(Total_Source, N, Properties)

        type(PropertiesType), intent(in) :: Properties
        type(NType), intent(in)          :: N

        real(kind = 8),    intent(inout) :: Total_Source

        integer :: i

        Total_Source = 0.0_8

        do i = 1,N%Element

            Total_Source = Total_Source + SUM(Properties%Elements(i)%Source)

        end do

    end subroutine Calculate_Spherical_Total_Source

    subroutine Calculate_Spherical_Scalar_Flux(Properties, N)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, k, m

        real(kind = 8), dimension(N%Angle) :: w

        call Calculate_w(w)

        do i = 1,N%Element

            Properties%Elements(i)%Scalar_Flux = 0.0_8

            do k = 1,N%Group

                do m = 1,N%Angle

                    Properties%Elements(i)%Scalar_Flux(k,:) = Properties%Elements(i)%Scalar_Flux(k,:) + 0.5_8*w(m)*Properties%Elements(i)%Flux(k,m,:)

                end do

            end do

        end do

    end subroutine Calculate_Spherical_Scalar_Flux

    subroutine Calculate_Spherical_Source(Properties, N, k_eff, i)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        real(kind = 8) :: k_eff

        integer, intent(in) :: i
        integer :: ang, g_index, g_s_index

        do ang = 1,N%Angle

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

    end subroutine Calculate_Spherical_Source

    subroutine Calculate_Half_Flux(Properties, N, Half_Flux, k_eff)

        type(PropertiesType), intent(inout)  :: Properties
        type(NType), intent(in)              :: N

        real(kind = 8), dimension(:,:,:), intent(inout) :: Half_Flux

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: T, K, M, F_out, F_in

        real(kind = 8), dimension(N%Degree+1) :: S, Source, Boundary_Flux

        integer :: i, g_index, g_s_index

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

        do g_index = 1, N%Group

            call Spherical_Boundary_Conditions(Properties, N, g_index, N%Angle, Boundary)

            Boundary_Flux(1) = Boundary

            do i = N%Element,1,-1

                dr = abs(Properties%Elements(i)%Coordinates(1,1) - Properties%Elements(i)%Coordinates(2,1))

                T = (-2.0_8/dr)*(-F_out - K) + Properties%Elements(i)%Sigma_t(g_index)*M

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

                    Source = Source + (2.0_8/dr)*MATMUL(F_in,Boundary_Flux)

                else

                    Source = Source + (2.0_8/dr)*MATMUL(F_in,Half_Flux(g_index,i+1,:))

                end if

                call Solve_Matrix(T, Source, Half_Flux(g_index,i,:))

            end do

        end do

    end subroutine Calculate_Half_Flux

    subroutine Construct_Spherical_Matrix(Properties, N, alpha, L, mu, w)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer :: i, g_index, ang

        real(kind = 8), dimension(N%Angle+1) :: alpha

        real(kind = 8), dimension(N%Element,N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: F_out

        real(kind = 8), dimension(:), intent(in) :: mu, w

        real(kind = 8) :: r_minus, r_plus, dr

        do i = 1, N%Element

            r_minus = Properties%Elements(i)%Coordinates(1,1)

            r_plus = Properties%Elements(i)%Coordinates(2,1)

            dr = abs(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))

            L(i,:,:) = Create_L(Properties,N,i)

            do g_index = 1, N%Group

                do ang = 1, N%Angle

                    if (ang > N%Angle/2) then

                        F_out = 0.0_8

                        F_out(1,1) = 1.0_8

                        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = (2.0_8*mu(ang)/dr)*(-F_out*r_minus**2 - Properties%Elements(i)%S_Matrix(1,:,:)) + (2.0_8/w(ang))*(2.0_8*alpha(ang)*L(i,:,:)) + Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix(:,:)

                    else

                        F_out = 0.0_8

                        F_out(2,2) = 1.0_8

                        Properties%Elements(i)%K_Matrix(g_index,ang,:,:) = (2.0_8*mu(ang)/dr)*(F_out*r_plus**2 - Properties%Elements(i)%S_Matrix(1,:,:)) + (2.0_8/w(ang))*(2.0_8*alpha(ang)*L(i,:,:)) + Properties%Elements(i)%Sigma_t(g_index)*Properties%Elements(i)%A_Matrix(:,:)

                    end if

                end do

            end do

        end do

    end subroutine Construct_Spherical_Matrix

    subroutine Up_Wind_Spherical_Source(Properties, N, g_index, ang, i, mu, w, alpha, L, Half_Flux)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                 :: i, ang, g_index

        real(kind = 8), dimension(:), intent(in) :: alpha

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: F_in

        real(kind = 8), intent(in) :: mu, w
        real(kind = 8), dimension(:,:), intent(in) :: L

        real(kind = 8), dimension(N%Degree+1) :: Half_Flux, Boundary_Flux

        real(kind = 8) :: r_minus, r_plus, dr, Boundary

        dr = abs(Properties%Elements(i)%Coordinates(1,1) - Properties%Elements(i)%Coordinates(2,1))

        F_in = 0.0_8

        call Spherical_Boundary_Conditions(Properties, N, g_index, ang, Boundary)

        Boundary_Flux = 0.0_8

        Boundary_Flux(1) = Boundary

        if (ang > N%Angle/2) then

            F_in(2,1) = 1.0_8

            r_plus = Properties%Elements(i)%Coordinates(2,1)

            if (i == N%Element) then

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(ang+1) + alpha(ang))*(2.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) - (2.0_8*mu/dr)*matmul(F_in,Boundary_Flux)*r_plus**2

            else

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(ang+1) + alpha(ang))*(2.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) - (2.0_8*mu/dr)*matmul(F_in,Properties%Elements(i+1)%Flux(g_index,ang,:))*r_plus**2

            end if

        else

            F_in(1,2) = 1.0_8

            r_minus = Properties%Elements(i)%Coordinates(1,1)

            if (i == 1) then

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(ang+1) + alpha(ang))*(2.0_8/w)*MATMUL(L(:,:),Half_Flux(:))

            else

                Properties%Elements(i)%Total_Source(g_index,ang,:) = Properties%Elements(i)%Source(g_index,ang,:) + (alpha(ang+1) + alpha(ang))*(2.0_8/w)*MATMUL(L(:,:),Half_Flux(:)) + (2.0_8*mu/dr)*matmul(F_in,Properties%Elements(i-1)%Flux(g_index,ang,:))*r_minus**2

            end if

        end if

    end subroutine Up_Wind_Spherical_Source

    subroutine Spherical_Boundary_Conditions(Properties, N, k, m, Right_Boundary)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in)                :: k, m

        real(kind = 8) :: Right_Boundary

        if (Properties%RBC == 1) then

            Right_Boundary = Properties%Elements(N%Element)%Flux(k,N%Angle+1-m,2)

        else if (Properties%RBC == 2) then

            Right_Boundary = 0.0_8

        else if (Properties%RBC == 3) then

            Right_Boundary = Properties%alpha*Properties%Elements(N%Element)%Flux(k,N%Angle+1-m,2)

        else if (Properties%RBC == 4) then

            print *, 'Periodic Boundary Conditions not possible'
            stop

        else if (Properties%RBC == 5) then

            Right_Boundary = Properties%Q_s

        end if

    end subroutine Spherical_Boundary_Conditions

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

                    M(index_1,index_2) = M(index_1,index_2) + w(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x_val)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x_val)!*Properties%Elements(i)%Det_Jacobian(j)

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

    Function Create_L(Properties,N,i) Result(L)

        type(PropertiesType) :: Properties
        type(NType)          :: N

        integer, intent(in) :: i
        integer             :: j, index_1, index_2

        real(kind = 8), dimension(N%Degree+1,N%Degree+1) :: L

        real(kind = 8), dimension(N%Degree+1) :: x, r
        real(kind = 8), dimension(N%Degree+1) :: w
        real(kind = 8)                        :: x_val

        L = 0.0_8

        call Generate_1D_Quad_Gauss_Points(N%Degree+1, x, w)

        r = (Properties%Elements(i)%Coordinates(1,1) + Properties%Elements(i)%Coordinates(2,1))/2.0_8 + x*(Properties%Elements(i)%Coordinates(2,1) - Properties%Elements(i)%Coordinates(1,1))/2.0_8

        do j = 1, N%Degree+1

            x_val = x(j)

            do index_1 = 1, N%Degree+1

                do index_2 = 1, N%Degree+1

                    L(index_1,index_2) = L(index_1,index_2) + w(j)*r(j)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_1, x_val)*Generate_Shape_Functions(Properties%Isoparametric_Coordinates(:,1), N, index_2, x_val)

                end do

            end do

        end do

    end Function Create_L

end module m_Spherical