module m_Sweep_Order_2D
!
! Purpose:
! To determine the order in which elements are swept in a 2D transport problem.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
    use m_Read_Properties
    use m_VTK_Reader
    use m_Calculate_mu_w
    use m_RZ

  implicit none

  contains

  subroutine determine_sweep_order_2D(Properties, N, Sweep_Order)

      type(PropertiesType), intent(inout) :: Properties
      type(NType), intent(in)             :: N
      integer, intent(out)                :: Sweep_Order(:,:)

      integer, dimension(:,:,:), allocatable :: Dependencies

      integer :: m, i, j, Sweep_Counter, iter, Total_Ordinates

      logical :: contains_negative, removed_element

      logical, dimension(N%Element) :: skip

      if(Properties%g == 0) allocate(Dependencies(N%Ordinates,N%Element,4))
      if(Properties%g == 1) allocate(Dependencies(N%Ordinates+N%Angle,N%Element,4))

      if(Properties%g == 0) Total_Ordinates = N%Ordinates
      if(Properties%g == 1) Total_Ordinates = N%Ordinates+N%Angle

      Sweep_Order = 0

      call Create_NEList(Properties, N, Dependencies)

      ! For each angle, determine the order in which the elements should be swept
      !$OMP PARALLEL DO PRIVATE(m, iter, i, j, Sweep_Counter, contains_negative, removed_element, skip)
      do m = 1, Total_Ordinates
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
                          exit
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
              print*, 'Error: Infinite loop in determine_sweep_order'
              stop
            end if
          end do
      end do
      !$OMP END PARALLEL DO
      
  end subroutine determine_sweep_order_2D

  subroutine Calculate_Unit_Vectors(Properties, N)

      type(PropertiesType), intent(inout) :: Properties
      type(NType), intent(in)             :: N

      integer :: i, j, node_num, k

      real(kind=8) :: dx, dy, mag

      do i = 1, N%Element
          allocate(Properties%Elements(i)%Unit_Vectors(Properties%Elements(i)%Number_of_Sides,2))
          do j = 1, Properties%Elements(i)%Number_of_Sides
              node_num = j + 2
              if (node_num > Properties%Elements(i)%Number_of_Sides) node_num = node_num - Properties%Elements(i)%Number_of_Sides
              if (node_num == Properties%Elements(i)%Number_of_Sides) then
                  dx = (Properties%Elements(i)%Coordinates(1,1) - Properties%Elements(i)%Coordinates(node_num,1))
                  dy = (Properties%Elements(i)%Coordinates(1,2) - Properties%Elements(i)%Coordinates(node_num,2))
                  Properties%Elements(i)%Side_Nodes(j,1) = node_num
                  Properties%Elements(i)%Side_Nodes(j,2) = 1
                  if (N%Degree > 1) then
                    Properties%Elements(i)%Side_Nodes(j,3) = Properties%Elements(i)%Number_of_Sides+1+(node_num-1)*(N%Degree-1)
                    do k = 3, N%Degree
                      Properties%Elements(i)%Side_Nodes(j,1+k) = Properties%Elements(i)%Side_Nodes(j,k) + 1
                    end do
                  end if
              else
                  dx = (Properties%Elements(i)%Coordinates(node_num+1,1) - Properties%Elements(i)%Coordinates(node_num,1))
                  dy = (Properties%Elements(i)%Coordinates(node_num+1,2) - Properties%Elements(i)%Coordinates(node_num,2))
                  Properties%Elements(i)%Side_Nodes(j,1) = node_num
                  Properties%Elements(i)%Side_Nodes(j,2) = node_num + 1
                  if (N%Degree > 1) then
                    Properties%Elements(i)%Side_Nodes(j,3) = Properties%Elements(i)%Number_of_Sides+1+(node_num-1)*(N%Degree-1)
                    do k = 3, N%Degree
                      Properties%Elements(i)%Side_Nodes(j,1+k) = Properties%Elements(i)%Side_Nodes(j,k) + 1
                    end do
                  end if
              end if
              mag = sqrt(dx**2 + dy**2)
              Properties%Elements(i)%Unit_Vectors(j,1) = dy/mag
              Properties%Elements(i)%Unit_Vectors(j,2) = -dx/mag
          end do
      end do

  end subroutine Calculate_Unit_Vectors

  ! Function to check if the nodes are in anticlockwise order if not order of nodes for dx, dy needs to be reversed
  logical function CheckAnticlockwise(Nodes, RefNode) result(anticlockwise)

      real(kind = 8), dimension(2,2)    :: Nodes
      real(kind = 8), dimension(2)      :: RefNode
      real(kind = 8)                    :: determinant, x1, y1, x2, y2, x0, y0

      x1 = Nodes(1,1)
      y1 = Nodes(1,2)
      x2 = Nodes(2,1)
      y2 = Nodes(2,2)
      x0 = RefNode(1)
      y0 = RefNode(2)

      determinant = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
      
      if (determinant > 0) then
          anticlockwise = .true.
      elseif (determinant < 0) then
          anticlockwise = .false.
      else
          print*, 'Nodes are collinear'
          stop
      end if

  end function CheckAnticlockwise

  subroutine Ordinates_2D(m,N_Angles,mu,eta,mu_val,eta_val)

      integer, intent(in)    :: m, N_Angles
      integer :: a, b
      integer :: i, j, counter, m_eff

      real(kind=8), intent(in), dimension(:) :: mu, eta

      real(kind=8) :: mu_val, eta_val

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

      if (m <= (N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = eta(b)
      else if (m >= (N_Angles+2)*N_Angles/8 + 1 .and. m <= 2*(N_Angles+2)*N_Angles/8) then
        mu_val = -mu(a)
        eta_val = eta(b)
      else if (m >= 2*(N_Angles+2)*N_Angles/8 + 1 .and.m <= 3*(N_Angles+2)*N_Angles/8) then
        mu_val = mu(a)
        eta_val = -eta(b)
      else if (m >= 3*(N_Angles+2)*N_Angles/8 + 1) then
        mu_val = -mu(a)
        eta_val = -eta(b)
      end if

  end subroutine Ordinates_2D

  subroutine Create_NEList(Properties, Nu, EEList)

    type(PropertiesType), intent(inout) :: Properties
    type(NType), intent(in) :: Nu

    integer :: i, j, n

    integer, dimension(Nu%N) :: node_counter

    integer, dimension(:,:), allocatable :: NEList

    integer, dimension(:,:,:), intent(inout) :: EEList

    allocate(NEList(Nu%N,8))

    node_counter = 0

    NEList = 0

    do i = 1, Nu%Element

      do j = 1, Properties%Elements(i)%Number_of_Sides

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

    integer :: i, j, k, l, m, ang, side_node_1, side_node_2, neighbour_side_node_1, neighbour_side_node_2, p, mu_counter

    integer, dimension(2) :: facet

    integer, dimension(:,:,:), intent(inout) :: EEList

    integer, dimension(:,:), allocatable :: EEList_1

    integer, dimension(Nu%Ordinates) :: p_indices, q_indices, w_indices, mu_indices

    real(kind=8), dimension(Nu%Angle/2) :: mu, eta

    real(kind = 8), dimension(Nu%Angle)  :: xsi

    real(kind=8) :: dot_product_value, mu_val, eta_val

    allocate(EEList_1(Nu%Element,4))

    EEList = 0

    EEList_1 = 0

    call Calculate_Unit_Vectors(Properties,Nu)

    call calculate_mu(mu)

    call calculate_mu(eta)

    mu = abs(mu)
    eta = abs(eta)

    do i = 1, Nu%Element

      do j = 1, Properties%Elements(i)%Number_of_Sides

        Properties%Elements(i)%Neighbours(j,1) = 0

        facet = 0

        do k = 1, 2

          facet(k) = Properties%Elements(i)%Cell_Pointers(mod(j+k,Properties%Elements(i)%Number_of_Sides)+1)

        end do

        do l = 1, size(NEList,2)

          do m = 1, size(NEList,2)

            if (NEList(facet(1),l) == NEList(facet(2),m) .and. NEList(facet(1),l) /= i .and. NEList(facet(1),l) /= 0) then

              EEList_1(i,j) = NEList(facet(1),l)

              Properties%Elements(i)%Neighbours(j,1) = EEList_1(i,j)

            end if

          end do

        end do

        if (EEList_1(i,j) /= 0) then

          do l = 1, Properties%Elements(EEList_1(i,j))%Number_of_Sides

            side_node_1 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,1))

            side_node_2 = Properties%Elements(i)%Cell_Pointers(Properties%Elements(i)%Side_Nodes(j,2))

            neighbour_side_node_1 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,1))

            neighbour_side_node_2 = Properties%Elements(EEList_1(i,j))%Cell_Pointers(Properties%Elements(EEList_1(i,j))%Side_Nodes(l,2))

            if (neighbour_side_node_1 == side_node_2 .and. neighbour_side_node_2 == side_node_1) then

              Properties%Elements(i)%Neighbours(j,2) = l

            end if

          end do

        else

          if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == 0.0_8 .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == 0.0_8) then
            Properties%Elements(i)%Neighbours(j,2) = 1
          else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),1) == Properties%Length(1) .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),1) == Properties%Length(1)) then
            Properties%Elements(i)%Neighbours(j,2) = 2
          else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == Properties%Length(2) .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == Properties%Length(2)) then
            Properties%Elements(i)%Neighbours(j,2) = 3
          else if (Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,1),2) == 0.0_8 .and. Properties%Elements(i)%Coordinates(Properties%Elements(i)%Side_Nodes(j,2),2) == 0.0_8) then
            Properties%Elements(i)%Neighbours(j,2) = 4
          end if

        end if

        allocate(Properties%Elements(i)%Sides(j)%Boundary(Nu%Group,Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes))

        if (Properties%g == 0) then
          if (Properties%Elements(i)%Neighbours(j,1) == 0) then
            allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))
          else
            allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes))
          end if
        else if (Properties%g == 1) then
          if (Properties%Elements(i)%Neighbours(j,1) == 0) then
            allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates+Nu%Angle,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(i)%Number_of_Nodes))
          else
            allocate(Properties%Elements(i)%Sides(j)%F_in_Matrix(Nu%Ordinates+Nu%Angle,Properties%Elements(i)%Number_of_Nodes,Properties%Elements(Properties%Elements(i)%Neighbours(j,1))%Number_of_Nodes))
          end if
        end if

      end do

    end do

    if (Properties%g == 0) then

    do ang = 1,Nu%Ordinates

      call Ordinates_2D(ang,Nu%Angle,mu,eta,mu_val,eta_val)

      do i = 1, Nu%Element

        do j = 1, Properties%Elements(i)%Number_of_Sides

          dot_product_value = Properties%Elements(i)%Unit_Vectors(j,1)*mu_val + Properties%Elements(i)%Unit_Vectors(j,2)*eta_val

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

    else if (Properties%g == 1) then

    call calculate_xsi(xsi)

    call Create_Ordinate_Arrays_RZ(Nu, p_indices, q_indices, w_indices, mu_indices)

    do ang = 1,Nu%Ordinates

      if (ang == 1) mu = -mu
      if (ang == Nu%Ordinates/4 + 1) mu = -mu
      if (ang == Nu%Ordinates/2 + 1) mu = -mu
      if (ang == 3*Nu%Ordinates/4 + 1) mu = -mu

      p = p_indices(ang)
      mu_counter = mu_indices(ang)

      do i = 1, Nu%Element

        do j = 1, Properties%Elements(i)%Number_of_Sides

          dot_product_value = Properties%Elements(i)%Unit_Vectors(j,1)*mu(mu_counter) + Properties%Elements(i)%Unit_Vectors(j,2)*xsi(p)

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

    do ang = 1,Nu%Angle

      do i = 1, Nu%Element

        do j = 1, Properties%Elements(i)%Number_of_Sides

          dot_product_value = -Properties%Elements(i)%Unit_Vectors(j,1)*SQRT(1.0_8 - xsi(ang)**2) + Properties%Elements(i)%Unit_Vectors(j,2)*xsi(ang)

          if (dot_product_value > 0) then

            EEList(Nu%Ordinates+ang,i,j) = EEList_1(i,j)

          elseif (dot_product_value == 0.0_8) then

            if(i > Properties%Elements(i)%Neighbours(j,1)) then
              EEList(Nu%Ordinates+ang,i,j) = EEList_1(i,j)
            else
              EEList(Nu%Ordinates+ang,i,j) = -EEList_1(i,j)
            end if

          elseif (dot_product_value < 0) then

            EEList(Nu%Ordinates+ang,i,j) = -EEList_1(i,j)

          end if
            
        end do

      end do

    end do

    end if

  end subroutine Create_EEList

end module m_Sweep_Order_2D
