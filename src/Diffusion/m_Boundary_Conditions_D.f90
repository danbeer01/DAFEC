module m_Boundary_Conditions_D
!
! Purpose:
! To apply boundary conditions to the matrix
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Gauss_Points
use m_Create_Hexahedral_Shape_Functions_D
use m_Create_Tetrahedral_Shape_Functions_D
use m_Create_Pyramidal_Shape_Functions_D
use m_Create_Prismatic_Shape_Functions_D
use m_Create_Quadrilateral_Shape_Functions_D
use m_Create_Triangular_Shape_Functions_D
use m_Create_Shape_Functions_D
use m_VTK_Reader

implicit none

contains

    subroutine Calculate_B_Matrix_1D(Properties,i,B_Matrix)

        type(PropertiesTypeD), intent(inout) :: Properties

        real(kind = 8), dimension(:,:), intent(inout) :: B_Matrix

        integer, intent(in) :: i

        integer :: Node_index

        integer :: k_left, k_right

        integer :: left_node, right_node

        logical :: left_boundary, right_boundary

        B_Matrix = 0.0_8

        k_left = 0
        k_right = 0

        left_node = 0
        right_node = 0

        left_boundary = .false.
        right_boundary = .false.

        do Node_index = 1, Properties%Elements(i)%Number_of_Nodes

            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - 0.0_8) < 1e-6) then
                k_left = k_left + 1
                left_node = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - Properties%Length(1)) < 1e-6) then
                k_right = k_right + 1
                right_node = Node_index
            end if

        end do

        if (k_left == 1) left_boundary = .true.
        if (k_right == 1) right_boundary = .true.

        if (left_boundary) then

            if (Properties%LBC == 2 .or. Properties%LBC == 3) then

                B_Matrix(left_node,left_node) = 0.5_8

                if (Properties%LBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%LBC == 0) then

                B_Matrix(left_node,left_node) = 1.0e9_8

            else if (Properties%LBC == 4) then

                B_Matrix(left_node,left_node) = 1.0e9_8

            end if

            if (Properties%g == 1 .or. Properties%g == 2) B_Matrix = 0.0_8

        end if

        if (right_boundary) then

            if (Properties%RBC == 2 .or. Properties%RBC == 3) then

                B_Matrix(right_node,right_node) = 0.5_8*(Properties%Length(1)**Properties%g)

                if (Properties%RBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%RBC == 0) then

                B_Matrix(right_node,right_node) = 1.0e9_8

            else if (Properties%RBC == 4) then

                B_Matrix(right_node,right_node) = 1.0e9_8

            end if

        end if

    end subroutine Calculate_B_Matrix_1D

    subroutine Calculate_B_Matrix_2D(Properties,N,i,B_Matrix)

        type(PropertiesTypeD), intent(inout) :: Properties
        type(NTypeD), intent(in)          :: N

        real(kind = 8), dimension(:,:), intent(inout) :: B_Matrix

        real(kind = 8) :: r

        integer, intent(in) :: i

        integer :: Node_index, j

        integer :: k_left, k_right, k_top, k_bottom

        integer, dimension((N%Degree+1)) :: left_nodes, right_nodes, top_nodes, bottom_nodes

        logical :: left_boundary, right_boundary, top_boundary, bottom_boundary

        B_Matrix = 0.0_8

        k_left = 0
        k_right = 0
        k_top = 0
        k_bottom = 0

        left_nodes = 0
        right_nodes = 0
        top_nodes = 0
        bottom_nodes = 0

        left_boundary = .false.
        right_boundary = .false.
        top_boundary = .false.
        bottom_boundary = .false.

        do Node_index = 1, Properties%Elements(i)%Number_of_Nodes

            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - 0.0_8) < 1e-6) then
                k_left = k_left + 1
                left_nodes(k_left) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - Properties%Length(1)) < 1e-6) then
                k_right = k_right + 1
                right_nodes(k_right) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,2) - 0.0_8) < 1e-6) then
                k_bottom = k_bottom + 1
                bottom_nodes(k_bottom) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,2) - Properties%Length(2)) < 1e-6) then
                k_top = k_top + 1
                top_nodes(k_top) = Node_index
            end if

        end do

        if (k_left == 1 + N%Degree) left_boundary = .true.
        if (k_right == 1 + N%Degree) right_boundary = .true.
        if (k_top == 1 + N%Degree) top_boundary = .true.
        if (k_bottom == 1 + N%Degree) bottom_boundary = .true.

        if (left_boundary) then

            if (Properties%LBC == 2 .or. Properties%LBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    B_Matrix = B_Matrix + Tri_Side_Boundary(Properties,N,i,left_nodes(1:k_left))

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
                    
                    B_Matrix = B_Matrix + Quad_Side_Boundary(Properties,N,i,left_nodes(1:k_left))

                end if

                if (Properties%LBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%LBC == 0) then

                do Node_index = 1, k_left

                    B_Matrix(left_nodes(Node_index),left_nodes(Node_index)) = 1.0e9_8

                end do

            end if

            if (Properties%g == 1) B_Matrix = 0.0_8

        end if

        if (right_boundary) then

            if (Properties%RBC == 2 .or. Properties%RBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    B_Matrix = B_Matrix + Tri_Side_Boundary(Properties,N,i,right_nodes(1:k_right))

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
                    
                    B_Matrix = B_Matrix + Quad_Side_Boundary(Properties,N,i,right_nodes(1:k_right))

                end if

                if (Properties%RBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%RBC == 0) then

                do Node_index = 1, k_right

                    B_Matrix(right_nodes(Node_index),right_nodes(Node_index)) = 1.0e9_8

                end do

            end if

            if (Properties%g == 1) B_Matrix = B_Matrix*Properties%Length(1)

        end if

        if (top_boundary) then

            if (Properties%TBC == 2 .or. Properties%TBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    B_Matrix = B_Matrix + Tri_Side_Boundary(Properties,N,i,top_nodes(1:k_top))

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
                    
                    B_Matrix = B_Matrix + Quad_Side_Boundary(Properties,N,i,top_nodes(1:k_top))

                end if

                if (Properties%TBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%TBC == 0) then

                do Node_index = 1, k_top

                    B_Matrix(top_nodes(Node_index),top_nodes(Node_index)) = 1.0e9_8

                end do

            end if

            r = 0.0_8

            if (Properties%g == 1) then
                do j = 1, k_top
                    r = r + Properties%Elements(i)%Coordinates(top_nodes(j),1)
                end do
                r = r/k_top
                B_Matrix = B_Matrix*r
            end if

        end if

        if (bottom_boundary) then

            if (Properties%BBC == 2 .or. Properties%BBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 5 .or. Properties%Elements(i)%Cell_Type == 22) then

                    B_Matrix = B_Matrix + Tri_Side_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                else if (Properties%Elements(i)%Cell_Type == 9 .or. Properties%Elements(i)%Cell_Type == 28) then
                    
                    B_Matrix = B_Matrix + Quad_Side_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                end if

                if (Properties%BBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%BBC == 0) then

                do Node_index = 1, k_bottom

                    B_Matrix(bottom_nodes(Node_index),bottom_nodes(Node_index)) = 1.0e9_8

                end do

            end if

            r = 0.0_8

            if (Properties%g == 1) then
                do j = 1, k_bottom
                    r = r + Properties%Elements(i)%Coordinates(bottom_nodes(j),1)
                end do
                r = r/k_bottom
                B_Matrix = B_Matrix*r
            end if

        end if

    end subroutine Calculate_B_Matrix_2D

    subroutine Calculate_B_Matrix_3D(Properties,N,i,B_Matrix)

        type(PropertiesTypeD), intent(inout) :: Properties
        type(NTypeD), intent(in)          :: N

        real(kind = 8), dimension(:,:), intent(inout) :: B_Matrix

        integer, intent(in) :: i

        integer :: Node_index

        integer :: k_left, k_right, k_top, k_bottom, k_front, k_back

        integer, dimension((N%Degree+1)**2) :: left_nodes, right_nodes, top_nodes, bottom_nodes, front_nodes, back_nodes

        logical :: left_boundary, right_boundary, top_boundary, bottom_boundary, front_boundary, back_boundary

        B_Matrix = 0.0_8

        k_left = 0
        k_right = 0
        k_top = 0
        k_bottom = 0
        k_front = 0
        k_back = 0

        left_nodes = 0
        right_nodes = 0
        top_nodes = 0
        bottom_nodes = 0
        front_nodes = 0
        back_nodes = 0

        left_boundary = .false.
        right_boundary = .false.
        top_boundary = .false.
        bottom_boundary = .false.
        front_boundary = .false.
        back_boundary = .false.

        do Node_index = 1, Properties%Elements(i)%Number_of_Nodes

            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - 0.0_8) < 1e-6) then
                k_left = k_left + 1
                left_nodes(k_left) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,1) - Properties%Length(1)) < 1e-6) then
                k_right = k_right + 1
                right_nodes(k_right) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,2) - 0.0_8) < 1e-6) then
                k_bottom = k_bottom + 1
                bottom_nodes(k_bottom) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,2) - Properties%Length(2)) < 1e-6) then
                k_top = k_top + 1
                top_nodes(k_top) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,3) - 0.0_8) < 1e-6) then
                k_front = k_front + 1
                front_nodes(k_front) = Node_index
            end if
            if (abs(Properties%Elements(i)%Coordinates(Node_index,3) - Properties%Length(3)) < 1e-6) then
                k_back = k_back + 1
                back_nodes(k_back) = Node_index
            end if

        end do

        if (k_left >= 2 + N%Degree) left_boundary = .true.
        if (k_right >= 2 + N%Degree) right_boundary = .true.
        if (k_top >= 2 + N%Degree) top_boundary = .true.
        if (k_bottom >= 2 + N%Degree) bottom_boundary = .true.
        if (k_front >= 2 + N%Degree) front_boundary = .true.
        if (k_back >= 2 + N%Degree) back_boundary = .true.

        if (left_boundary) then

            if (Properties%LBC == 2 .or. Properties%LBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,left_nodes(1:k_left))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,left_nodes(1:k_left))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,left_nodes(1:k_left))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,left_nodes(1:k_left))

                end if

                if (Properties%LBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%LBC == 0) then

                do Node_index = 1, k_left

                    B_Matrix(left_nodes(Node_index),left_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

        if (right_boundary) then

            if (Properties%RBC == 2 .or. Properties%RBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,right_nodes(1:k_right))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,right_nodes(1:k_right))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,right_nodes(1:k_right))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,right_nodes(1:k_right))

                end if

                if (Properties%RBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%RBC == 0) then

                do Node_index = 1, k_right

                    B_Matrix(right_nodes(Node_index),right_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

        if (top_boundary) then

            if (Properties%TBC == 2 .or. Properties%TBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,top_nodes(1:k_top))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,top_nodes(1:k_top))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,top_nodes(1:k_top))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,top_nodes(1:k_top))

                end if

                if (Properties%TBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%TBC == 0) then

                do Node_index = 1, k_top

                    B_Matrix(top_nodes(Node_index),top_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

        if (bottom_boundary) then

            if (Properties%BBC == 2 .or. Properties%BBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,bottom_nodes(1:k_bottom))

                end if

                if (Properties%BBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%BBC == 0) then

                do Node_index = 1, k_bottom

                    B_Matrix(bottom_nodes(Node_index),bottom_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

        if (front_boundary) then

            if (Properties%FBC == 2 .or. Properties%FBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,front_nodes(1:k_front))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,front_nodes(1:k_front))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,front_nodes(1:k_front))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,front_nodes(1:k_front))

                end if

                if (Properties%FBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%FBC == 0) then

                do Node_index = 1, k_front

                    B_Matrix(front_nodes(Node_index),front_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

        if (back_boundary) then

            if (Properties%BaBC == 2 .or. Properties%BaBC == 3) then

                if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

                    B_Matrix = B_Matrix + Hex_Face_Boundary(Properties,N,i,back_nodes(1:k_back))

                else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
                    
                    B_Matrix = B_Matrix + Tet_Face_Boundary(Properties,N,i,back_nodes(1:k_back))

                else if (Properties%Elements(i)%Cell_Type == 13) then

                    B_Matrix = B_Matrix + Pris_Face_Boundary(Properties,N,i,back_nodes(1:k_back))

                else if (Properties%Elements(i)%Cell_Type == 14) then
    
                    B_Matrix = B_Matrix + Pyr_Face_Boundary(Properties,N,i,back_nodes(1:k_back))

                end if

                if (Properties%BaBC == 3) B_Matrix = ((1-Properties%alpha)/(1+Properties%alpha))*B_Matrix

            else if (Properties%BaBC == 0) then

                do Node_index = 1, k_back

                    B_Matrix(back_nodes(Node_index),back_nodes(Node_index)) = 1.0e9_8

                end do

            end if

        end if

    end subroutine Calculate_B_Matrix_3D

    subroutine Periodic_Boundary_3D(N,Properties,Periodic_Pairs)
            
        type(NTypeD), intent(in)          :: N
        type(PropertiesTypeD), intent(in) :: Properties

        integer             :: i, j, k, counter, counter_1, counter_2, counter_3, counter_4, counter_5, counter_6

        integer, dimension(:,:), allocatable          :: Periodic_Nodes
        real(kind = 8), dimension(:,:,:), allocatable :: Periodic_Nodes_Coordinates
        integer, dimension(:,:), allocatable          :: Periodic_Pairs

        allocate(Periodic_Nodes(6,NINT((real(N%N))**(2.0/3.0))))
        allocate(Periodic_Nodes_Coordinates(6,NINT((real(N%N))**(2.0/3.0)),2))

        Periodic_Nodes = 0

        counter_1 = 0
        counter_2 = 0
        counter_3 = 0
        counter_4 = 0
        counter_5 = 0
        counter_6 = 0

        do i = 1, N%Element
            do j = 1, size(Properties%Elements(i)%Coordinates,1)
                if (Properties%LBC == 4 .and. Properties%RBC == 4) then
                    if (Properties%Elements(i)%Coordinates(j,1) == 0.0_8) then
                        if (any(Periodic_Nodes(1,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_1 = counter_1 + 1
                            Periodic_Nodes(1,counter_1) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(1,counter_1,:) = Properties%Elements(i)%Coordinates(j,2:3)
                        end if
                    else if (Properties%Elements(i)%Coordinates(j,1) == Properties%Length(1)) then
                        if (any(Periodic_Nodes(2,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_2 = counter_2 + 1
                            Periodic_Nodes(2,counter_2) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(2,counter_2,:) = Properties%Elements(i)%Coordinates(j,2:3)
                        end if
                    end if
                else if (Properties%LBC == 4 .or. Properties%RBC == 4) then
                    write(*,*) 'Error: Periodic Boundary Conditions must be applied to both left and right boundaries'
                    stop
                end if
                if (Properties%BBC == 4 .and. Properties%TBC == 4) then
                    if (Properties%Elements(i)%Coordinates(j,2) == 0.0_8) then
                        if (any(Periodic_Nodes(3,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_3 = counter_3 + 1
                            Periodic_Nodes(3,counter_3) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(3,counter_3,1) = Properties%Elements(i)%Coordinates(j,1)
                            Periodic_Nodes_Coordinates(3,counter_3,2) = Properties%Elements(i)%Coordinates(j,3)
                        end if
                    else if (Properties%Elements(i)%Coordinates(j,2) == Properties%Length(2)) then
                        if (any(Periodic_Nodes(4,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_4 = counter_4 + 1
                            Periodic_Nodes(4,counter_4) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(4,counter_4,1) = Properties%Elements(i)%Coordinates(j,1)
                            Periodic_Nodes_Coordinates(4,counter_4,2) = Properties%Elements(i)%Coordinates(j,3)
                        end if
                    end if
                else if (Properties%TBC == 4 .or. Properties%BBC == 4) then
                    write(*,*) 'Error: Periodic Boundary Conditions must be applied to both top and bottom boundaries'
                    stop
                end if
                if (Properties%FBC == 4 .and. Properties%BaBC == 4) then
                    if (Properties%Elements(i)%Coordinates(j,3) == 0.0_8) then
                        if (any(Periodic_Nodes(5,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_5 = counter_5 + 1
                            Periodic_Nodes(5,counter_5) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(5,counter_5,:) = Properties%Elements(i)%Coordinates(j,1:2)
                        end if
                    else if (Properties%Elements(i)%Coordinates(j,3) == Properties%Length(3)) then
                        if (any(Periodic_Nodes(6,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_6 = counter_6 + 1
                            Periodic_Nodes(6,counter_6) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(6,counter_6,:) = Properties%Elements(i)%Coordinates(j,1:2)
                        end if
                    end if
                else if (Properties%BaBC == 4 .or. Properties%FBC == 4) then
                    write(*,*) 'Error: Periodic Boundary Conditions must be applied to both front and back boundaries'
                    stop
                end if
            end do
        end do

        if (counter_1 /= counter_2) then
            write(*,*) 'Error: Periodic Boundary Conditions not possible due to different number of nodes on the left and right boundaries'
            stop
        end if
        if (counter_3 /= counter_4) then
            write(*,*) 'Error: Periodic Boundary Conditions not possible due to different number of nodes on the top and bottom boundaries'
            stop
        end if
        if (counter_5 /= counter_6) then
            write(*,*) 'Error: Periodic Boundary Conditions not possible due to different number of nodes on the front and back boundaries'
            stop
        end if

        allocate(Periodic_Pairs(counter_1+counter_3+counter_5,2))

        counter = 0

        if (Properties%LBC == 4 .and. Properties%RBC == 4) then
            do j = 1, counter_1
                do k = 1, counter_1
                    if (ABS(Periodic_Nodes_Coordinates(1,j,1) - Periodic_Nodes_Coordinates(2,k,1)) < 1e-5) then
                        if (ABS(Periodic_Nodes_Coordinates(1,j,2) - Periodic_Nodes_Coordinates(2,k,2)) < 1e-5) then
                            counter = counter + 1
                            Periodic_Pairs(counter,:) = [Periodic_Nodes(1,j),Periodic_Nodes(2,k)]
                        end if
                    end if
                end do
            end do
        end if
        if (Properties%BBC == 4 .and. Properties%TBC == 4) then
            do j = 1, counter_3
                do k = 1, counter_3
                    if (ABS(Periodic_Nodes_Coordinates(3,j,1) - Periodic_Nodes_Coordinates(4,k,1)) < 1e-5) then
                        if (ABS(Periodic_Nodes_Coordinates(3,j,2) - Periodic_Nodes_Coordinates(4,k,2)) < 1e-5) then
                            counter = counter + 1
                            Periodic_Pairs(counter,:) = [Periodic_Nodes(3,j),Periodic_Nodes(4,k)]
                        end if
                    end if
                end do
            end do
        end if
        if (Properties%FBC == 4 .and. Properties%BaBC == 4) then
            do j = 1, counter_5
                do k = 1, counter_5
                    if (ABS(Periodic_Nodes_Coordinates(5,j,1) - Periodic_Nodes_Coordinates(6,k,1)) < 1e-5) then
                        if (ABS(Periodic_Nodes_Coordinates(5,j,2) - Periodic_Nodes_Coordinates(6,k,2)) < 1e-5) then
                            counter = counter + 1
                            Periodic_Pairs(counter,:) = [Periodic_Nodes(5,j),Periodic_Nodes(6,k)]
                        end if
                    end if
                end do
            end do
        end if

        Periodic_Pairs = Periodic_Pairs - 1

    end subroutine Periodic_Boundary_3D

    subroutine Periodic_Boundary_2D(N,Properties,Periodic_Pairs)
            
        type(NTypeD), intent(in)          :: N
        type(PropertiesTypeD), intent(in) :: Properties

        integer             :: i, j, k, counter, counter_1, counter_2, counter_3, counter_4

        integer, dimension(:,:), allocatable        :: Periodic_Nodes
        real(kind = 8), dimension(:,:), allocatable :: Periodic_Nodes_Coordinates
        integer, dimension(:,:), allocatable        :: Periodic_Pairs

        allocate(Periodic_Nodes(4,NINT(sqrt(real(N%N)))))
        allocate(Periodic_Nodes_Coordinates(4,NINT(sqrt(real(N%N)))))

        Periodic_Nodes = 0

        counter_1 = 0
        counter_2 = 0
        counter_3 = 0
        counter_4 = 0

        do i = 1, N%Element
            do j = 1, size(Properties%Elements(i)%Coordinates,1)
                if (Properties%LBC == 4 .and. Properties%RBC == 4) then
                    if (Properties%Elements(i)%Coordinates(j,1) == 0.0_8) then
                        if (any(Periodic_Nodes(1,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_1 = counter_1 + 1
                            Periodic_Nodes(1,counter_1) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(1,counter_1) = Properties%Elements(i)%Coordinates(j,2)
                        end if
                    else if (Properties%Elements(i)%Coordinates(j,1) == Properties%Length(1)) then
                        if (any(Periodic_Nodes(2,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_2 = counter_2 + 1
                            Periodic_Nodes(2,counter_2) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(2,counter_2) = Properties%Elements(i)%Coordinates(j,2)
                        end if
                    end if
                else if (Properties%LBC == 4 .or. Properties%RBC == 4) then
                    write(*,*) 'Error: Periodic Boundary Conditions must be applied to both left and right boundaries'
                    stop
                end if
                if (Properties%BBC == 4 .and. Properties%TBC == 4) then
                    if (Properties%Elements(i)%Coordinates(j,2) == 0.0_8) then
                        if (any(Periodic_Nodes(3,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_3 = counter_3 + 1
                            Periodic_Nodes(3,counter_3) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(3,counter_3) = Properties%Elements(i)%Coordinates(j,1)
                        end if
                    else if (Properties%Elements(i)%Coordinates(j,2) == Properties%Length(2)) then
                        if (any(Periodic_Nodes(4,:) == Properties%Elements(i)%Cell_Pointers(j))) then
                            continue
                        else
                            counter_4 = counter_4 + 1
                            Periodic_Nodes(4,counter_4) = Properties%Elements(i)%Cell_Pointers(j)
                            Periodic_Nodes_Coordinates(4,counter_4) = Properties%Elements(i)%Coordinates(j,1)
                        end if
                    end if
                else if (Properties%TBC == 4 .or. Properties%BBC == 4) then
                    write(*,*) 'Error: Periodic Boundary Conditions must be applied to both top and bottom boundaries'
                    stop
                end if
            end do
        end do

        if (counter_1 /= counter_2) then
            write(*,*) 'Error: Periodic Boundary Conditions not possible due to different number of nodes on the left and right boundaries'
            stop
        end if
        if (counter_3 /= counter_4) then
            write(*,*) 'Error: Periodic Boundary Conditions not possible due to different number of nodes on the top and bottom boundaries'
            stop
        end if

        allocate(Periodic_Pairs(counter_1+counter_3,2))

        counter = 0

        if (Properties%LBC == 4 .and. Properties%RBC == 4) then
            do j = 1, counter_1
                do k = 1, counter_1
                    if (ABS(Periodic_Nodes_Coordinates(1,j) - Periodic_Nodes_Coordinates(2,k)) < 1e-5) then
                        counter = counter + 1
                        Periodic_Pairs(counter,:) = [Periodic_Nodes(1,j),Periodic_Nodes(2,k)]
                    end if
                end do
            end do
        end if
        if (Properties%BBC == 4 .and. Properties%TBC == 4) then
            do j = 1, counter_3
                do k = 1, counter_3
                    if (ABS(Periodic_Nodes_Coordinates(3,j) - Periodic_Nodes_Coordinates(4,k)) < 1e-5) then
                        counter = counter + 1
                        Periodic_Pairs(counter,:) = [Periodic_Nodes(3,j),Periodic_Nodes(4,k)]
                    end if
                end do
            end do
        end if

        Periodic_Pairs = Periodic_Pairs - 1

    end subroutine Periodic_Boundary_2D

end module