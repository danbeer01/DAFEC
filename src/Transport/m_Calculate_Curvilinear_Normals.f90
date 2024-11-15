module m_Calculate_Curvilinear_Normals
use m_Read_Properties
use m_Gauss_Points
use m_Create_Quadrilateral_Shape_Functions
implicit none

contains

    subroutine Calculate_Curvilinear_Unit_Vectors(Properties, N, i, j, unit_vector)

        type(PropertiesType), intent(inout) :: Properties
        type(NType), intent(in)             :: N

        integer, intent(in) :: i, j

        real(kind=8), dimension(:), allocatable :: xi, w

        integer :: Num_Gauss_Points, k, l, p

        integer, dimension(:), allocatable :: Nodes

        real(kind = 8), dimension(:), allocatable   :: dx_dxi, dy_dxi

        real(kind = 8), dimension(:,:), allocatable :: normal_vector

        real(kind = 8), dimension(:,:) :: unit_vector

        real(kind = 8) :: magnitude

        Num_Gauss_Points = 2*N%Degree+2

        allocate(dx_dxi(Num_Gauss_Points), dy_dxi(Num_Gauss_Points))

        allocate(xi(Num_Gauss_Points), w(Num_Gauss_Points))

        allocate(normal_vector(Num_Gauss_Points, 2))

        allocate(Nodes(N%Degree+1))

        call Generate_1D_Quad_Gauss_Points(Num_Gauss_Points, xi, w)

        Nodes(1) = 1
        Nodes(2) = 2
        do l = 3, N%Degree+1
            if (l == 3) then
                Nodes(l) = 5
            else
                Nodes(l) = 1 + Nodes(l-1)
            end if
        end do

        dx_dxi = 0.0_8
        dy_dxi = 0.0_8

        do l = 1, Num_Gauss_Points

            do p = 1, size(Nodes)
                k = Properties%Elements(i)%Side_Nodes(j,p)
                dx_dxi(l) = dx_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,1)
                dy_dxi(l) = dy_dxi(l) + Generate_Shape_Functions_Derivative_xi(Properties, N%Degree, Nodes(p), 1.0_8, xi(l))*Properties%Elements(i)%Coordinates(k,2)
            end do

            normal_vector(l,:) = [-dy_dxi(l), dx_dxi(l)]

            magnitude = sqrt(normal_vector(l,1)**2 + normal_vector(l,2)**2)

            normal_vector(l,:) = normal_vector(l,:)/magnitude

            unit_vector(l,:) = normal_vector(l,:)

        end do

        ! mean_unit_vector = sum(normal_vector, dim=1)

        ! magnitude = sqrt(mean_unit_vector(1)**2 + mean_unit_vector(2)**2)

        ! mean_unit_vector = mean_unit_vector/magnitude

    end subroutine Calculate_Curvilinear_Unit_Vectors

end module m_Calculate_Curvilinear_Normals