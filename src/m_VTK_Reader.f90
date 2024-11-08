module m_VTK_Reader
!
! Purpose:
! To read in a VTK file and store the data in a mesh object
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
implicit none

type Mesh
    integer :: N_nodes
    integer :: N_cells
    integer :: N_degree
    integer :: Max_Nodes_Per_Element
    integer :: Dimension
    real(8), allocatable :: Nodes(:,:)
    integer, allocatable :: Cell_Type(:)
    integer, allocatable :: Cell_Pointers(:,:)
    integer, allocatable :: N_Cell_Pointers(:)
    integer, allocatable :: Material_ID(:)
contains
    PROCEDURE,pass :: read_VTK_file
    PROCEDURE,pass :: Fix_Cell_Pointers
    PROCEDURE,pass :: Determine_Degree
end type Mesh

contains

    subroutine read_VTK_file(this)

        class(Mesh), intent(out) :: this

        character(100)    :: line
        integer           :: i, j, Max_Nodes_Per_Element, Nodes_Per_Element

        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*) line, this%N_nodes

        allocate(this%nodes(this%N_nodes,3))

        do i = 1, this%N_nodes
            read(1,*) this%nodes(i,:)
        end do

        read(1,*)
        read(1,*) line, this%N_cells

        Max_Nodes_Per_Element = 0

        do i = 1, this%N_cells
            read(1,*) Nodes_Per_Element
            if (Nodes_Per_Element > Max_Nodes_Per_Element) Max_Nodes_Per_Element = Nodes_Per_Element
        end do

        this%Max_Nodes_Per_Element = Max_Nodes_Per_Element

        allocate(this%cell_pointers(this%N_cells,Max_Nodes_Per_Element))
        allocate(this%N_cell_pointers(this%N_cells))

        rewind(1)

        do i = 1, 7 + this%N_nodes
            read(1,*)
        end do

        do i = 1, this%N_cells
            read(1, *) this%N_cell_pointers(i), (this%cell_pointers(i,j), j = 1, this%N_cell_pointers(i))
        end do

        this%cell_pointers = this%cell_pointers + 1

        read(1,*)
        read(1,*)

        allocate(this%cell_type(this%N_cells))

        do i = 1, this%N_cells
            read(1,*) this%cell_type(i)
            if (i == 1) call this%Determine_Degree()
            if(this%N_degree == 3 .and. this%cell_type(i) == 12) call Fix_Cell_Pointers(this,i)
        end do

        read(1,*)
        read(1,*)
        read(1,*)
        read(1,*)

        allocate(this%material_ID(this%N_cells))

        do i = 1, this%N_cells
            read(1,*) this%material_ID(i)
        end do

        close(1)

    end subroutine read_VTK_file

    subroutine Fix_Cell_Pointers(this,element)

        class(Mesh), intent(inout) :: this

        integer, dimension(64) :: Old_Cell_Pointers, New_Cell_Pointers

        integer :: i

        integer :: element

        Old_Cell_Pointers = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 14, 24, 25, 10, 11, 26, 27, 12, 13, 28, 29, 30, 31, 16, 17, 23, 22, 18, 19, 20, 21, 32, 35, 34, 33, 40, 41, 42, 43, 53, 54, 55, 52, 44, 45, 46, 47, 48, 49, 50, 51, 36, 37, 38, 39, 56, 57, 58, 59, 60, 61, 62, 63]

        Old_Cell_Pointers = Old_Cell_Pointers + 1

        do i = 1, 64
            New_Cell_Pointers(Old_Cell_Pointers(i)) = this%cell_pointers(element,i)
        end do

        this%cell_pointers(element,:) = New_Cell_Pointers

    end subroutine Fix_Cell_Pointers

    subroutine Determine_Degree(this)

        class(Mesh), intent(inout) :: this

        if (this%Cell_Type(1) == 3 .or. this%Cell_Type(1) == 21) then

            this%Dimension = 1

        else if (this%Cell_Type(1) == 5 .or. this%Cell_Type(1) == 22 .or. this%Cell_Type(1) == 9 .or. this%Cell_Type(1) == 28) then

            this%Dimension = 2

        else if (this%Cell_Type(1) == 12 .or. this%Cell_Type(1) == 29 .or. this%Cell_Type(1) == 10 .or. this%Cell_Type(1) == 24 .or. this%Cell_Type(1) == 13 .or. this%Cell_Type(1) == 14) then

            this%Dimension = 3

        else

            write(*,*) 'Error: Cannot determine dimension of mesh'
            stop

        end if

        if (this%Dimension == 1) then

            this%N_degree = this%Max_Nodes_Per_Element - 1

        else if (this%Dimension == 2) then

            if (this%Max_Nodes_Per_Element /= 36 .and. mod(0.5*(-1.0 + sqrt(1.0 + 8.0*real(this%Max_Nodes_Per_Element))), 1.0) == 0.0) then

                this%N_degree = INT(0.5*(-1.0 + sqrt(1.0 + 8.0*real(this%Max_Nodes_Per_Element)))) - 1

            else if (mod(sqrt(real(this%Max_Nodes_Per_Element)), 1.0) == 0.0) then

                this%N_degree = INT(sqrt(real(this%Max_Nodes_Per_Element))) - 1

            end if

        else if (this%Dimension == 3) then

            if (mod((real(this%Max_Nodes_Per_Element))**(1.0/3.0), 1.0) == 0.0) then

                this%N_degree = NINT((real(this%Max_Nodes_Per_Element))**(1.0/3.0)) - 1

            else if (mod((NINT(10000.0*((3.0*real(this%Max_Nodes_Per_Element) + SQRT(9.0*real(this%Max_Nodes_Per_Element)**2 - 1.0/27.0))**(1.0/3.0) +  (3.0*real(this%Max_Nodes_Per_Element) - SQRT(9.0*real(this%Max_Nodes_Per_Element)**2 - 1.0/27.0))**(1.0/3.0)))), 10000) == 0.0) then

                this%N_degree = NINT((3.0*real(this%Max_Nodes_Per_Element) + SQRT(9.0*real(this%Max_Nodes_Per_Element)**2 - 1.0/27.0))**(1.0/3.0) +  (3.0*real(this%Max_Nodes_Per_Element) - SQRT(9.0*real(this%Max_Nodes_Per_Element)**2 - 1.0/27.0))**(1.0/3.0)) - 2

            else if (this%Max_Nodes_Per_Element == 15 .or. this%Max_Nodes_Per_Element == 14) then

                this%N_degree = 2

            else

                write(*,*) 'Error: Cannot determine degree of polynomial, assumed linear'
                this%N_degree = 1

            end if

        end if


    end subroutine Determine_Degree

end module m_VTK_Reader
