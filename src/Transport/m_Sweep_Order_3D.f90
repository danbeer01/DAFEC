module m_Sweep_Order_3D
!
! Purpose:
! To determine the order in which elements are swept in a 3D transport problem.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
  use m_Read_Properties
  use m_VTK_Reader
  use m_Calculate_mu_w

  implicit none
  
  contains
  
  subroutine determine_sweep_order_3D(Properties, N, Sweep_Order)

      type(PropertiesType), intent(inout) :: Properties
      type(NType), intent(in)             :: N
      integer, intent(out)                :: Sweep_Order(:,:)

      integer, dimension(:,:,:), allocatable :: Dependencies

      integer :: m, i, j, Sweep_Counter, iter

      logical :: contains_negative, removed_element

      logical, dimension(N%Element) :: skip

      allocate(Dependencies(N%Ordinates,N%Element,6))

      Sweep_Order = 0

      call Create_NEList(Properties, N, Dependencies)

      ! For each angle, determine the order in which the elements should be swept
      !$OMP PARALLEL DO PRIVATE(m, iter, i, j, Sweep_Counter, contains_negative, removed_element, skip)
      do m = 1, N%Ordinates
          ! Reset the sweep order
          Sweep_Counter = 0
          skip = .false.
          iter = 0
          ! Iterate until all elements have been added to the sweep order
          do
            iter = iter + 1
            do i = 1, N%Element
                ! Check if the element has already been added to the sweep order
                if (.not. skip(i)) then
                  ! Assume the element has all positive dependencies
                  contains_negative = .false.
                  ! Go through every dependency of the element
                  do j = 1, Properties%Elements(i)%Number_of_Sides
                      ! If any negative dependencies are found, check if they have already been added to the sweep order
                      if (Dependencies(m,i,j) < 0) then
                        ! If sweep order is empty then the negative dependency is not in the sweep order
                        if (Sweep_Counter == 0) then
                          contains_negative = .true.
                          exit
                        end if
                        removed_element = .false.
                        ! Go through the sweep order to check if the negative dependency has already been added
                        if (skip(-Dependencies(m,i,j))) then
                          removed_element = .true.
                          Dependencies(m,i,j) = 0
                        end if
                        ! If the negative dependency has not been added to the sweep order then the element has a negative dependency
                        if (.not. removed_element) then
                          contains_negative = .true.
                          ! exit
                        end if
                      end if
                  end do
                  ! If the element has no negative dependencies then add it to the sweep order
                  if (.not. contains_negative) then
                    Sweep_Counter = Sweep_Counter + 1
                    Sweep_Order(m,Sweep_Counter) = i
                    ! Mark the element as added to the sweep order so it will be skipped for any future iterations
                    skip(i) = .true.
                  end if
                end if
            end do
            ! If all elements have been added to the sweep order then exit the loop
            if (Sweep_Counter == N%Element) exit
            if (iter > N%Element) then
              print*, 'Error: Infinite loop in determine_sweep_order', Sweep_Counter
              stop
            end if
          end do
      end do
      !$OMP END PARALLEL DO
      
  end subroutine determine_sweep_order_3D

  subroutine Calculate_Unit_Vectors(Properties, N)

      type(PropertiesType), intent(inout) :: Properties
      type(NType), intent(in)             :: N

      integer :: i, j

      real(kind=8), dimension(3) :: vector1, vector2, normal

      real(kind=8) :: mag

      do i = 1, N%Element
          allocate(Properties%Elements(i)%Unit_Vectors(Properties%Elements(i)%Number_of_Sides,3))
          do j = 1, Properties%Elements(i)%Number_of_Sides
            if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then
              if (j == 1) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,4,9,10,11,12,25/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,4,9,10,11,12,13,14,15,16,33,34,35,36/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 2) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/5,6,7,8/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/5,6,7,8,13,14,15,16,26/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/5,6,7,8,17,18,19,20,21,22,23,24,37,38,39,40/)
                vector1 = Properties%Elements(i)%Coordinates(6,:) - Properties%Elements(i)%Coordinates(5,:)
                vector2 = Properties%Elements(i)%Coordinates(7,:) - Properties%Elements(i)%Coordinates(5,:)
              else if (j == 3) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,6,5/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,6,5,9,18,13,17,23/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,6,5,9,10,27,28,18,17,26,25,41,42,43,44/)
                vector1 = Properties%Elements(i)%Coordinates(6,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 4) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/4,3,7,8/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/4,3,7,8,11,19,15,20,24/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/4,3,7,8,14,13,29,30,21,22,32,31,50,49,52,51/)
                vector1 = Properties%Elements(i)%Coordinates(8,:) - Properties%Elements(i)%Coordinates(4,:)
                vector2 = Properties%Elements(i)%Coordinates(7,:) - Properties%Elements(i)%Coordinates(4,:)
              else if (j == 5) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,5,8,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,5,8,4,17,16,20,12,21/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/1,5,8,4,25,26,24,23,32,31,15,16,54,55,56,53/)
                vector1 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(8,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 6) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/2,6,7,3/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/2,6,7,3,18,14,19,10,22/)
                if(N%Degree == 3) Properties%Elements(i)%Side_Nodes(j,:) = (/2,6,7,3,27,28,19,20,30,29,12,11,45,48,47,46/)
                vector1 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(6,:)
                vector2 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(6,:)
              end if
            else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then
              if (j == 1) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,5,6,7/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 2) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,4,5,9,8/)
                vector1 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 3) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,4,6,10,9/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(2,:)
                vector2 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(2,:)
              else if (j == 4) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,3,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,3,4,7,10,8/)
                vector1 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
              end if
            else if (Properties%Elements(i)%Cell_Type == 13) then
              if (j == 1) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,6,5/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,6,5,8,12,14,11/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(2,:)
                vector2 = Properties%Elements(i)%Coordinates(6,:) - Properties%Elements(i)%Coordinates(2,:)
              else if (j == 2) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/3,1,4,6/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/3,1,4,6,9,10,15,12/)
                vector1 = Properties%Elements(i)%Coordinates(6,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 3) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,5,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,5,4,7,11,13,10/)
                vector1 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 4) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,7,8,9,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 5) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/4,5,6,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/4,5,6,13,14,15,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(4,:)
                vector2 = Properties%Elements(i)%Coordinates(6,:) - Properties%Elements(i)%Coordinates(4,:)
              end if
            else if (Properties%Elements(i)%Cell_Type == 14) then
              if (j == 1) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,4/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,3,4,6,7,8,9,14/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 2) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,5,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,2,5,6,11,10,0,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(2,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(1,:)
              else if (j == 3) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,5,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/2,3,5,7,12,11,0,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(3,:) - Properties%Elements(i)%Coordinates(2,:)
                vector2 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(2,:)
              else if (j == 4) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/3,4,5,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/3,4,5,8,13,12,0,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(3,:)
                vector2 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(3,:)
              else if (j == 5) then
                if(N%Degree == 1) Properties%Elements(i)%Side_Nodes(j,:) = (/1,4,5,0/)
                if(N%Degree == 2) Properties%Elements(i)%Side_Nodes(j,:) = (/1,4,5,9,13,10,0,0,0/)
                vector1 = Properties%Elements(i)%Coordinates(5,:) - Properties%Elements(i)%Coordinates(1,:)
                vector2 = Properties%Elements(i)%Coordinates(4,:) - Properties%Elements(i)%Coordinates(1,:)
              end if
            end if
            normal(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
            normal(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3)
            normal(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)
            mag = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
            Properties%Elements(i)%Unit_Vectors(j,:) = normal/mag
          end do
      end do

  end subroutine Calculate_Unit_Vectors

  ! Function to check if the nodes are in anticlockwise order if not order of nodes for dx, dy needs to be reversed
  logical function CheckAnticlockwise(Nodes) result(anticlockwise)

      real(kind = 8), dimension(5,3)    :: Nodes
      real(kind = 8)                    :: x1, y1, z1, x2, y2, z2, x0, y0, z0
      real(kind = 8), dimension(3)      :: vector1, vector2, normal

      x1 = Nodes(1,1)
      y1 = Nodes(1,2)
      z1 = Nodes(1,3)
      x2 = Nodes(2,1)
      y2 = Nodes(2,2)
      z2 = Nodes(2,3)

      x0 = Nodes(3,1)
      y0 = Nodes(3,2)
      z0 = Nodes(3,3)

      vector1 = (/x1 - x0, y1 - y0, z1 - z0/)
      vector2 = (/x2 - x0, y2 - y0, z2 - z0/)

      normal(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
      normal(2) = vector1(3)*vector2(1) - vector1(1)*vector2(3)
      normal(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)

      if (sum(normal) > 0) then
          anticlockwise = .false.
      elseif (sum(normal) < 0) then
          anticlockwise = .true.
      else
          print*, 'Nodes are collinear'
          stop
      end if

  end function CheckAnticlockwise

  subroutine Ordinates(m,N_Angles,mu,mu_val,eta_val,xi_val)

      integer, intent(in)    :: m, N_Angles
      integer :: a, b, c
      integer :: i, j, counter, m_eff

      real(kind=8), intent(in), dimension(:) :: mu

      real(kind=8) :: mu_val, eta_val, xi_val

      m_eff = MOD(m,(N_Angles+2)*N_Angles/8)

      if (m_eff == 0) m_eff = (N_Angles+2)*N_Angles/8

      counter = 0
      do j = 1, N_Angles/2
        do i = 1, N_Angles/2
          a = i
          b = j
          counter = counter + 1
          if (counter == m_eff) exit
          if (i == N_Angles/2-j+1) exit
        end do
        if (counter == m_eff) exit
      end do

      c = N_Angles/2 + 2 - a - b

      if (m <= (N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = mu(b)
        xi_val = mu(c)
      else if (m >= (N_Angles+2)*N_Angles/8 + 1 .and. m <= 2*(N_Angles+2)*N_Angles/8) then
        mu_val = -mu(a)
        eta_val = mu(b)
        xi_val = mu(c)
      else if (m >= 2*(N_Angles+2)*N_Angles/8 + 1 .and.m <= 3*(N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = -mu(b)
        xi_val = mu(c)
      else if (m >= 3*(N_Angles+2)*N_Angles/8 + 1 .and. m <= 4*(N_Angles+2)*N_Angles/8) then
        mu_val = -mu(a)
        eta_val = -mu(b)
        xi_val = mu(c)
      else if (m >= 4*(N_Angles+2)*N_Angles/8 + 1 .and. m <= 5*(N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = mu(b)
        xi_val = -mu(c)
      else if (m >= 5*(N_Angles+2)*N_Angles/8 + 1 .and. m <= 6*(N_Angles+2)*N_Angles/8) then
        mu_val = -mu(a)
        eta_val = mu(b)
        xi_val = -mu(c)
      else if (m >= 6*(N_Angles+2)*N_Angles/8 + 1 .and. m <= 7*(N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = -mu(b)
        xi_val = -mu(c)
      else if (m >= 7*(N_Angles+2)*N_Angles/8 + 1 .and. m <= (N_Angles+2)*N_Angles) then
        mu_val = -mu(a)
        eta_val = -mu(b)
        xi_val = -mu(c)
      end if

  end subroutine Ordinates

  subroutine Create_NEList(Properties, Nu, EEList)

    type(PropertiesType), intent(inout) :: Properties
    type(NType), intent(in) :: Nu

    integer :: i, j, n

    integer, dimension(Nu%N) :: node_counter

    integer, dimension(:,:), allocatable :: NEList

    integer, dimension(:,:,:), intent(inout) :: EEList

    allocate(NEList(Nu%N,50))

    node_counter = 0

    NEList = 0

    do i = 1, Nu%Element

      do j = 1, Properties%Elements(i)%Number_of_Nodes

        n = Properties%Elements(i)%Cell_Pointers(j)

        node_counter(n) = node_counter(n) + 1

        NEList(n,node_counter(n)) = i

      end do

    end do

    call Create_EEList(Properties,Nu,NEList,EEList)

  end subroutine Create_NEList

  subroutine Create_EEList(Properties,Nu,NEList,EEList)

    type(PropertiesType), intent(inout) :: Properties
    type(NType), intent(in) :: Nu

    integer, dimension(:,:), intent(in) :: NEList

    integer :: i, j, k, l, ang

    integer :: i_1, i_2, i_3, i_4

    integer :: side_node_1, side_node_2, side_node_3, side_node_4

    integer :: neighbour_side_node_1, neighbour_side_node_2, neighbour_side_node_3, neighbour_side_node_4

    integer, dimension(4) :: facet

    integer, dimension(:,:,:), intent(inout) :: EEList

    integer, dimension(:,:), allocatable :: EEList_1

    real(kind=8), dimension(Nu%Angle/2) :: mu

    real(kind=8) :: dot_product_value, mu_val, eta_val, xi_val

    allocate(EEList_1(Nu%Element,6))

    EEList = 0

    EEList_1 = 0

    call Calculate_Unit_Vectors(Properties,Nu)

    call calculate_mu(mu)

    do i = 1, Nu%Element

      Properties%Elements(i)%Neighbours = 0

      do j = 1, Properties%Elements(i)%Number_of_Sides

        Properties%Elements(i)%Neighbours(j,1) = 0

        facet = 0

        do k = 1, Properties%Elements(i)%Sides(j)%Nodes_on_Side

          facet(k) = Properties%Elements(i)%Cell_Pointers(Create_facet(Properties,i,j,k))

        end do

        if (Properties%Elements(i)%Sides(j)%Nodes_on_Side == 4) then

          if(Nu%Degree == 1) allocate(Properties%Elements(i)%Sides(j)%Neighbour_Nodes(4))
          if(Nu%Degree == 2) allocate(Properties%Elements(i)%Sides(j)%Neighbour_Nodes(8))

          do i_1 = 1, size(NEList,2)

            if (NEList(facet(1),i_1) /= i .and. NEList(facet(1),i_1) /= 0) then

              do i_2 = 1, size(NEList,2)

                do i_3 = 1, size(NEList,2)

                  do i_4 = 1, size(NEList,2)

                    if (NEList(facet(1),i_1) == NEList(facet(2),i_2) .and. NEList(facet(2),i_2) == NEList(facet(3),i_3) .and. NEList(facet(3),i_3) == NEList(facet(4),i_4)) then

                      EEList_1(i,j) = NEList(facet(1),i_1)

                      Properties%Elements(i)%Neighbours(j,1) = EEList_1(i,j)

                    end if

                  end do

                end do

              end do

            end if

          end do

          if (EEList_1(i,j) /= 0) then

            do l = 1, Properties%Elements(EEList_1(i,j))%Number_of_Sides

              if (Properties%Elements(EEList_1(i,j))%Sides(l)%Nodes_on_Side == 4) then
  
                side_node_1 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,1))
    
                side_node_2 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,2))
    
                side_node_3 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,3))
    
                side_node_4 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,4))
    
                neighbour_side_node_1 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,1))
    
                neighbour_side_node_2 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,2))
    
                neighbour_side_node_3 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,3))
    
                neighbour_side_node_4 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,4))
    
                if (side_node_1 == neighbour_side_node_1 .and. side_node_2 == neighbour_side_node_2 .and. side_node_3 == neighbour_side_node_3 .and. side_node_4 == neighbour_side_node_4) then
    
                  Properties%Elements(i)%Neighbours(j,2) = l

                else if (side_node_1 == neighbour_side_node_2 .and. side_node_2 == neighbour_side_node_1 .and. side_node_3 == neighbour_side_node_4 .and. side_node_4 == neighbour_side_node_3) then

                  Properties%Elements(i)%Neighbours(j,2) = l

                  if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,1,4,5/)
                  if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,1,4,5,7,10,13,11/)
    
                end if

              end if
  
            end do
 
          else
  
            if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == Properties%Left_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == Properties%Left_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),1) == Properties%Left_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),1) == Properties%Left_B) then
              Properties%Elements(i)%Neighbours(j,2) = 1
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == Properties%Right_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == Properties%Right_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),1) == Properties%Right_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),1) == Properties%Right_B) then
              Properties%Elements(i)%Neighbours(j,2) = 2
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == Properties%Top_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == Properties%Top_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),2) == Properties%Top_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),2) == Properties%Top_B) then
              Properties%Elements(i)%Neighbours(j,2) = 3
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == Properties%Bottom_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == Properties%Bottom_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),2) == Properties%Bottom_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),2) == Properties%Bottom_B) then
              Properties%Elements(i)%Neighbours(j,2) = 4
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),3) == Properties%Back_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),3) == Properties%Back_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),3) == Properties%Back_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),3) == Properties%Back_B) then
              Properties%Elements(i)%Neighbours(j,2) = 5
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),3) == Properties%Front_B.and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),3) == Properties%Front_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),3) == Properties%Front_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,4),3) == Properties%Front_B) then
              Properties%Elements(i)%Neighbours(j,2) = 6
            end if
  
          end if

        else if (Properties%Elements(i)%Sides(j)%Nodes_on_Side == 3) then

          if(Nu%Degree == 1) allocate(Properties%Elements(i)%Sides(j)%Neighbour_Nodes(3))
          if(Nu%Degree == 2) allocate(Properties%Elements(i)%Sides(j)%Neighbour_Nodes(6))

          do i_1 = 1, size(NEList,2)

            if (NEList(facet(1),i_1) /= i .and. NEList(facet(1),i_1) /= 0) then

              do i_2 = 1, size(NEList,2)

                do i_3 = 1, size(NEList,2)

                  if (NEList(facet(1),i_1) == NEList(facet(2),i_2) .and. NEList(facet(2),i_2) == NEList(facet(3),i_3)) then

                    EEList_1(i,j) = NEList(facet(1),i_1)

                    Properties%Elements(i)%Neighbours(j,1) = EEList_1(i,j)

                  end if

                end do

              end do

            end if

          end do

          if (EEList_1(i,j) /= 0) then

            do l = 1, Properties%Elements(EEList_1(i,j))%Number_of_Sides
  
              side_node_1 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,1))
  
              side_node_2 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,2))
  
              side_node_3 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,3))
   
              neighbour_side_node_1 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,1))
  
              neighbour_side_node_2 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,2))
  
              neighbour_side_node_3 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,3))
    
              if (side_node_1 == neighbour_side_node_1 .and. side_node_2 == neighbour_side_node_2 .and. side_node_3 == neighbour_side_node_3) then
  
                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/1,2,3/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/1,2,3,5,6,7/)

              else if (side_node_1 == neighbour_side_node_2 .and. side_node_2 == neighbour_side_node_3 .and. side_node_3 == neighbour_side_node_1) then
  
                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/3,1,2/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/3,1,2,7,5,6/)

              else if (side_node_1 == neighbour_side_node_3 .and. side_node_2 == neighbour_side_node_1 .and. side_node_3 == neighbour_side_node_2) then
  
                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,3,1/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,3,1,6,7,5/)

              else if (side_node_1 == neighbour_side_node_2 .and. side_node_2 == neighbour_side_node_1 .and. side_node_3 == neighbour_side_node_3) then

                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,1,3/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/2,1,3,5,7,6/)

              else if (side_node_1 == neighbour_side_node_3 .and. side_node_2 == neighbour_side_node_2 .and. side_node_3 == neighbour_side_node_1) then

                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/3,2,1/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/3,2,1,6,5,7/)

              else if (side_node_1 == neighbour_side_node_1 .and. side_node_2 == neighbour_side_node_3 .and. side_node_3 == neighbour_side_node_2) then

                Properties%Elements(i)%Neighbours(j,2) = l

                if(Nu%Degree == 1) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/1,3,2/)
                if(Nu%Degree == 2) Properties%Elements(i)%Sides(j)%Neighbour_Nodes = (/1,3,2,7,6,5/)
  
              end if
  
            end do
  
          else
  
            if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == Properties%Left_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == Properties%Left_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),1) == Properties%Left_B) then
              Properties%Elements(i)%Neighbours(j,2) = 1
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == Properties%Right_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == Properties%Right_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),1) == Properties%Right_B) then
              Properties%Elements(i)%Neighbours(j,2) = 2
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == Properties%Top_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == Properties%Top_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),2) == Properties%Top_B) then
              Properties%Elements(i)%Neighbours(j,2) = 3
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == Properties%Bottom_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == Properties%Bottom_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),2) == Properties%Bottom_B) then
              Properties%Elements(i)%Neighbours(j,2) = 4
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),3) == Properties%Back_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),3) == Properties%Back_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),3) == Properties%Back_B) then
              Properties%Elements(i)%Neighbours(j,2) = 5
            else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),3) == Properties%Front_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),3) == Properties%Front_B .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,3),3) == Properties%Front_B) then
              Properties%Elements(i)%Neighbours(j,2) = 6
            end if
  
          end if

        end if

        allocate(Properties%Elements(i)%Sides(j)%Boundary(Nu%Group,Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes))

        if (Properties%Elements(i)%Neighbours(j,1) == 0) then
          allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))
        else
          allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes))
        end if

      end do

    end do

    do ang = 1,Nu%Ordinates

      call Ordinates(ang,Nu%Angle,mu,mu_val,eta_val,xi_val)

      do i = 1, Nu%Element

        do j = 1, Properties%Elements(i)%Number_of_Sides

          dot_product_value = Properties%Elements(i)%Unit_Vectors(j,1)*mu_val + Properties%Elements(i)%Unit_Vectors(j,2)*eta_val + Properties%Elements(i)%Unit_Vectors(j,3)*xi_val

          if(abs(dot_product_value) < 1e-10) dot_product_value = 0.0_8

          if (dot_product_value > 0) then

            EEList(ang,i,j) = EEList_1(i,j)

          elseif (dot_product_value == 0.0_8) then

            if(i > Properties%Elements(i)%Neighbours(j,1)) then
              EEList(ang,i,j) = EEList_1(i,j)
            else
              EEList(ang,i,j) = -EEList_1(i,j)
            end if

          elseif (dot_product_value < 0) then

            EEList(ang,i,j) = -EEList_1(i,j)

          end if
            
        end do

      end do

    end do

  end subroutine Create_EEList

  function Create_facet(Properties,i,j,k) result (node)

      type(PropertiesType), intent(in) :: Properties
      
      integer, intent(in) :: i, j, k
    
      integer :: node

      integer, dimension(:), allocatable :: nodes

      allocate(nodes(Properties%Elements(i)%Sides(j)%Nodes_on_Side))

      if (Properties%Elements(i)%Cell_Type == 12 .or. Properties%Elements(i)%Cell_Type == 29) then

        if (j == 1) then
          nodes = (/1,2,3,4/)
        else if (j == 2) then
          nodes = (/5,6,7,8/)
        else if (j == 3) then
          nodes = (/1,2,6,5/)
        else if (j == 4) then
          nodes = (/4,3,7,8/)
        else if (j == 5) then
          nodes = (/1,5,8,4/)
        else if (j == 6) then
          nodes = (/2,6,7,3/)
        end if

      else if (Properties%Elements(i)%Cell_Type == 10 .or. Properties%Elements(i)%Cell_Type == 24) then

        if (j == 1) then
          nodes = (/1,2,3/)
        else if (j == 2) then
          nodes = (/1,2,4/)
        else if (j == 3) then
          nodes = (/2,3,4/)
        else if (j == 4) then
          nodes = (/1,3,4/)
        end if

      else if (Properties%Elements(i)%Cell_Type == 13) then

        if (j == 1) then
          nodes = (/2,3,6,5/)
        else if (j == 2) then
          nodes = (/3,1,4,6/)
        else if (j == 3) then
          nodes = (/1,2,5,4/)
        else if (j == 4) then
          nodes = (/1,2,3,0/)
        else if (j == 5) then
          nodes = (/4,5,6,0/)
        end if

      else if (Properties%Elements(i)%Cell_Type == 14) then

        if (j == 1) then
          nodes = (/1,2,3,4/)
        else if (j == 2) then
          nodes = (/1,2,5,0/)
        else if (j == 3) then
          nodes = (/2,3,5,0/)
        else if (j == 4) then
          nodes = (/3,4,5,0/)
        else if (j == 5) then
          nodes = (/1,4,5,0/)
        end if

      end if

      node = nodes(k)
  
    end function Create_facet

end module m_Sweep_Order_3D
