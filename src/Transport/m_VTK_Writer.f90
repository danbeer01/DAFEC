module m_VTK_Writer
!
! Purpose:
! To create a vtk file for visualisation of the neutron flux from my data
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties
use m_Results
use m_Constants_Mod
use m_VTK_Reader

contains

    subroutine Output_VTK(this, Results, N)
        ! Subroutine to create .vtk file
        type(Mesh), intent(inout)         :: this
        class(ResultsType), intent(inout) :: Results
        class(NType), intent(in)          :: N

        integer                              :: ii
        integer                              :: jj

        open(6, file = 'vtk_out.vtk', status = 'replace', action = 'write')

        write(6,1) "# vtk DataFile Version 2.0"
        write(6,2) "outputfile"
        write(6,3) "ASCII"
        write(6,4) "DATASET UNSTRUCTURED_GRID"
        write(6,5) "POINTS", this%N_nodes, "double"

        do ii = 1, this%N_nodes
            write(6,*) this%nodes(ii,1), this%nodes(ii,2), this%nodes(ii,3)
        end do

        ! Insert mesh
        write(6,*)

        this%Cell_Pointers = this%Cell_Pointers - 1

        write(6,7) "CELLS", this%N_cells, SUM(this%N_Cell_Pointers) + this%N_cells
        do ii = 1, this%N_cells
            write(6,8) this%N_Cell_Pointers(ii), this%Cell_Pointers(ii,1:this%N_Cell_Pointers(ii)) 
        end do

        ! Insert cell type
        write(6,*)
        write(6,9) "CELL_TYPES", this%N_cells
        do ii = 1, this%N_cells
            write(6,10) this%Cell_Type(ii)
        end do

        ! Insert cell data
        write(6,*)
        write(6,11) "CELL_DATA", this%N_cells
        write(6,*)
        write(6,12) "SCALARS blocks int 1"
        write(6,12) "LOOKUP_TABLE default"
        do ii = 1, this%N_cells
            write(6,13) ii-1
        end do

        ! Insert materials
        write(6,*)
        write(6, 17) "SCALARS material int 1"
        write(6,12) "LOOKUP_TABLE default"
        do ii = 1, this%N_cells
            write(6,13) this%Material_ID(ii)
        end do
        
        ! Insert neutron flux for each group
        do jj = 1, N%Group
            write(6,*)
            if (jj == 1) write(6,16) "POINT_DATA", this%N_nodes 
            write(6,14) "SCALARS scalar_flux_", jj, " double 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, this%N_nodes
                if (Results%Scalar_Flux(jj,ii) .lt. 1.0E-9_dp) then
                    Results%Scalar_Flux(jj,ii) = 0.0_dp
                end if
                write(6,15) Results%Scalar_Flux(jj,ii)
            end do
        end do

        do kk = 1, N%Ordinates
            do jj = 1, N%Group
                write(6,*)
                write(6,14) "SCALARS angular_flux_", kk, " double 1"
                write(6,12) "LOOKUP_TABLE default"
                do ii = 1, this%N_nodes
                    ! if (Results%Angular_Flux(jj,ii,kk) .lt. 1.0E-9_dp) then
                    !     Results%Angular_Flux(jj,ii,kk) = 0.0_dp
                    ! end if
                    write(6,*) Results%Angular_Flux(jj,ii,kk)
                end do
            end do
        end do

        1 format(A26)
        2 format(A10)
        3 format(A5)
        4 format(A25)
        5 format(A6, 1X, I0, 1X, A6)
        ! 6 format (ES20.14, 1X, ES20.14, 1X, ES20.14)
        7 format (A5, 1X, I0, 1X, I0)
        8 format(16(I0, 1X), I0)
        9 format (A10, 1X, I0)
        10 format(I2)
        11 format(A9, 1X, I0)
        12 format(A20)
        13 format(I0)
        14 format(A20, I0, A9)
        15 format(ES20.14)
        16 format(A10, 1X, I0)
        17 format(A22)
        ! 18 format(A21)
        ! 19 format(A31, I0, A9)
        ! 20 format(A34, I0, A9)
        close(6)

    end subroutine Output_VTK

    subroutine Output_Higher_Order_VTK(Properties, this, Results, N)
        ! Subroutine to create .vtk file
        type(Mesh), intent(inout)         :: this
        class(ResultsType), intent(inout) :: Results
        class(NType), intent(in)          :: N
        type(PropertiesType), intent(in)  :: Properties

        integer              :: ii
        integer              :: jj
        integer              :: counter
        integer              :: New_Number_of_Cells
        integer, allocatable :: New_Cell_Pointers(:,:)
        integer, allocatable :: Node_Numberings(:,:)

        open(6, file = 'vtk_out.vtk', status = 'replace', action = 'write')

        write(6,1) "# vtk DataFile Version 2.0"
        write(6,2) "outputfile"
        write(6,3) "ASCII"
        write(6,4) "DATASET UNSTRUCTURED_GRID"
        write(6,5) "POINTS", this%N_nodes, "double"

        do ii = 1, this%N_nodes
            write(6,6) this%nodes(ii,1), this%nodes(ii,2), this%nodes(ii,3)
        end do

        ! Insert mesh
        write(6,*)

        this%Cell_Pointers = this%Cell_Pointers - 1

        if (this%Cell_Type(1) == 9) then

            New_Number_of_Cells = this%N_cells*(N%Degree**2)

            allocate(New_Cell_Pointers(New_Number_of_Cells, 4))

            allocate(Node_Numberings(N%Degree**2, 4))

            call Calculate_Node_Numberings(Properties, N%Degree, Node_Numberings)

            counter = 1
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    New_Cell_Pointers(counter, 1) = this%Cell_Pointers(ii, Node_Numberings(jj, 1))
                    New_Cell_Pointers(counter, 2) = this%Cell_Pointers(ii, Node_Numberings(jj, 2))
                    New_Cell_Pointers(counter, 3) = this%Cell_Pointers(ii, Node_Numberings(jj, 3))
                    New_Cell_Pointers(counter, 4) = this%Cell_Pointers(ii, Node_Numberings(jj, 4))
                    counter = counter + 1
                end do
            end do

            write(6,7) "CELLS", New_Number_of_Cells, New_Number_of_Cells*5
            do ii = 1, New_Number_of_Cells
                write(6,8) 4, New_Cell_Pointers(ii,:)
            end do

            ! Insert cell type
            write(6,*)
            write(6,9) "CELL_TYPES", New_Number_of_Cells
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    write(6,10) 9
                end do
            end do

            ! Insert cell data
            write(6,*)
            write(6,11) "CELL_DATA", New_Number_of_Cells
            write(6,*)
            write(6,12) "SCALARS blocks int 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, New_Number_of_Cells
                write(6,13) ii-1
            end do

            ! Insert materials
            write(6,*)
            write(6, 17) "SCALARS material int 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    write(6,13) this%Material_ID(ii)
                end do
            end do

        else if (this%Cell_Type(1) == 5) then

            if (N%Degree > 3) then
                write(*,*) "ERROR: Only up to third order triangular meshes are supported"
                stop
            end if

            New_Number_of_Cells = this%N_cells*(N%Degree**2)

            allocate(New_Cell_Pointers(New_Number_of_Cells, 3))

            allocate(Node_Numberings(N%Degree**2, 3))

            Node_Numberings(:,1) = (/1,4,5,9,10,8,10,6,7/)
            Node_Numberings(:,2) = (/4,5,2,10,6,7,9,10,8/)
            Node_Numberings(:,3) = (/9,10,6,8,7,3,4,5,10/)

            counter = 1
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    New_Cell_Pointers(counter, 1) = this%Cell_Pointers(ii, Node_Numberings(jj, 1))
                    New_Cell_Pointers(counter, 2) = this%Cell_Pointers(ii, Node_Numberings(jj, 2))
                    New_Cell_Pointers(counter, 3) = this%Cell_Pointers(ii, Node_Numberings(jj, 3))
                    counter = counter + 1
                end do
            end do

            write(6,7) "CELLS", New_Number_of_Cells, New_Number_of_Cells*4
            do ii = 1, New_Number_of_Cells
                write(6,8) 3, New_Cell_Pointers(ii,:)
            end do

            ! Insert cell type
            write(6,*)
            write(6,9) "CELL_TYPES", New_Number_of_Cells
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    write(6,10) 5
                end do
            end do

            ! Insert cell data
            write(6,*)
            write(6,11) "CELL_DATA", New_Number_of_Cells
            write(6,*)
            write(6,12) "SCALARS blocks int 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, New_Number_of_Cells
                write(6,13) ii-1
            end do

            ! Insert materials
            write(6,*)
            write(6, 17) "SCALARS material int 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, this%N_cells
                do jj = 1, N%Degree**2
                    write(6,13) this%Material_ID(ii)
                end do
            end do

        else

            close(6)

            write(*,*) "ERROR: Only 2D triangular and quadrilateral meshes are supported"
            stop

        end if
        
        ! Insert neutron flux for each group
        do kk = 1, (N%Angle+2)*N%Angle/2
            do jj = 1, N%Group
                write(6,*)
                if (kk == 1) write(6,16) "POINT_DATA", this%N_nodes 
                write(6,14) "SCALARS angular_flux_", kk, " double 1"
                write(6,12) "LOOKUP_TABLE default"
                do ii = 1, this%N_nodes
                    if (Results%Angular_Flux(jj,ii,kk) .lt. 1.0E-9_dp) then
                        Results%Angular_Flux(jj,ii,kk) = 0.0_dp
                    end if
                    write(6,15) Results%Angular_Flux(jj,ii,kk)
                end do
            end do
        end do

        do jj = 1, N%Group
            write(6,*)
            ! if (kk == 1) write(6,16) "POINT_DATA", this%N_nodes 
            write(6,14) "SCALARS scalar_flux_", jj, " double 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, this%N_nodes
                if (Results%Scalar_Flux(jj,ii) .lt. 1.0E-9_dp) then
                    Results%Scalar_Flux(jj,ii) = 0.0_dp
                end if
                write(6,15) Results%Scalar_Flux(jj,ii)
            end do
        end do


        1 format(A26)
        2 format(A10)
        3 format(A5)
        4 format(A25)
        5 format(A6, 1X, I0, 1X, A6)
        6 format (ES20.14, 1X, ES20.14, 1X, ES20.14)
        7 format (A5, 1X, I0, 1X, I0)
        8 format(16(I0, 1X), I0)
        9 format (A10, 1X, I0)
        10 format(I2)
        11 format(A9, 1X, I0)
        12 format(A20)
        13 format(I0)
        14 format(A20, I0, A9)
        15 format(ES20.14)
        16 format(A10, 1X, I0)
        17 format(A22)
        ! 18 format(A21)
        ! 19 format(A31, I0, A9)
        ! 20 format(A34, I0, A9)
        close(6)

    end subroutine Output_Higher_Order_VTK

    subroutine Output_Discontinuous_VTK(Properties, this, N)
        ! Subroutine to create .vtk file
        type(PropertiesType), intent(in)  :: Properties
        type(Mesh), intent(inout)         :: this
        class(NType), intent(in)          :: N

        integer                              :: ii
        integer                              :: jj
        integer                              :: kk
        integer                              :: ll
        integer                              :: Counter
        integer                              :: Number_of_Nodes
        integer, allocatable                 :: New_Cell_Pointers(:,:)

        allocate(New_Cell_Pointers(N%Element,(N%Degree+1)**3))

        Number_of_Nodes = sum([(Properties%Elements(ii)%Number_of_Nodes, ii = 1, N%Element)])

        open(6, file = 'vtk_out.vtk', status = 'replace', action = 'write')

        write(6,1) "# vtk DataFile Version 2.0"
        write(6,2) "outputfile"
        write(6,3) "ASCII"
        write(6,4) "DATASET UNSTRUCTURED_GRID"
        write(6,5) "POINTS", Number_of_Nodes, "double"

        Counter = 0
        do ii = 1, N%Element
            do jj = 1, Properties%Elements(ii)%Number_of_Nodes
                write(6,*) Properties%Elements(ii)%Coordinates(jj,1), Properties%Elements(ii)%Coordinates(jj,2), Properties%Elements(ii)%Coordinates(jj,3)
                New_Cell_Pointers(ii,jj) = Counter
                Counter = Counter + 1
            end do
        end do

        ! Insert mesh
        write(6,*)

        write(6,7) "CELLS", this%N_cells, SUM(this%N_Cell_Pointers) + this%N_cells
        do ii = 1, this%N_cells
            write(6,8) this%N_Cell_Pointers(ii), New_Cell_Pointers(ii,1:this%N_Cell_Pointers(ii)) 
        end do

        ! Insert cell type
        write(6,*)
        write(6,9) "CELL_TYPES", this%N_cells
        do ii = 1, this%N_cells
            write(6,10) this%Cell_Type(ii)
        end do

        ! Insert cell data
        write(6,*)
        write(6,11) "CELL_DATA", this%N_cells
        write(6,*)
        write(6,12) "SCALARS blocks int 1"
        write(6,12) "LOOKUP_TABLE default"
        do ii = 1, this%N_cells
            write(6,13) ii-1
        end do

        ! Insert materials
        write(6,*)
        write(6, 17) "SCALARS material int 1"
        write(6,12) "LOOKUP_TABLE default"
        do ii = 1, this%N_cells
            write(6,13) this%Material_ID(ii)
        end do
        
        ! Insert neutron flux for each group
        do kk = 1, N%Ordinates
            do jj = 1, N%Group
                write(6,*)
                if (kk == 1) write(6,16) "POINT_DATA", Number_of_Nodes 
                write(6,14) "SCALARS angular_flux_", kk, " double 1"
                write(6,12) "LOOKUP_TABLE default"
                do ii = 1, N%Element
                    do ll = 1,Properties%Elements(ii)%Number_of_Nodes
                        ! if (Properties%Elements(ii)%Flux(jj,kk,ll) .lt. 1.0E-9_dp) then
                        !     write(6,15) 0.0_dp
                        ! else
                        !     write(6,15) Properties%Elements(ii)%Flux(jj,kk,ll)
                        ! end if
                        write(6,*) Properties%Elements(ii)%Flux(jj,kk,ll)
                    end do
                end do
            end do
        end do

        do jj = 1, N%Group
            write(6,*)
            write(6,14) "SCALARS scalar_flux_", jj, " double 1"
            write(6,12) "LOOKUP_TABLE default"
            do ii = 1, N%Element
                do ll = 1,Properties%Elements(ii)%Number_of_Nodes
                    if (Properties%Elements(ii)%Scalar_Flux(jj,ll) .lt. 1.0E-9_dp) then
                        write(6,15) 0.0_dp
                    else
                        write(6,15) Properties%Elements(ii)%Scalar_Flux(jj,ll)
                    end if
                end do
            end do
        end do

        1 format(A26)
        2 format(A10)
        3 format(A5)
        4 format(A25)
        5 format(A6, 1X, I0, 1X, A6)
        ! 6 format (ES20.14, 1X, ES20.14, 1X, ES20.14)
        7 format (A5, 1X, I0, 1X, I0)
        8 format(16(I0, 1X), I0)
        9 format (A10, 1X, I0)
        10 format(I2)
        11 format(A9, 1X, I0)
        12 format(A20)
        13 format(I0)
        14 format(A20, I0, A9)
        15 format(ES20.14)
        16 format(A10, 1X, I0)
        17 format(A22)
        ! 18 format(A21)
        ! 19 format(A31, I0, A9)
        ! 20 format(A34, I0, A9)
        close(6)

    end subroutine Output_Discontinuous_VTK

    subroutine Calculate_Node_Numberings(Properties, Degree, Node_Numberings)

        type(PropertiesType), intent(in) :: Properties
        integer, intent(in)              :: Degree

        real(kind = 8) :: Normal_Coordinates((Degree+1)**2, 2)
        real(kind = 8) :: x, y

        integer :: Node_Numberings(:,:), Conversion_Vector((Degree+1)**2)
        integer :: i, j, Coord_index, Counter

        Coord_index = 1

        do i = 0,Degree

            do j = 0,Degree

                x = 1.0_8 - 2.0_8*j/real(Degree,kind=8)

                y = 1.0_8 - 2.0_8*i/real(Degree,kind=8)

                Normal_Coordinates(Coord_index,:) = [x, y]

                Coord_index = Coord_index + 1

            end do

        end do

        do i = 1,(Degree+1)**2

            do j = 1,(Degree+1)**2

                if (ABS(Properties%Isoparametric_Coordinates(j,1) - Normal_Coordinates(i,1)) < 1e-5 .and. ABS(Properties%Isoparametric_Coordinates(j,2) - Normal_Coordinates(i,2)) < 1e-5) then

                    Conversion_Vector(i) = j

                end if

            end do

        end do

        Counter = 1

        do i = 1,Degree**2

            Node_Numberings(i,1) = Conversion_Vector(Counter)
            if (mod((real(Counter)+1.0)/(Real(Degree)+1.0), 1.0) == 0) Counter = Counter + 1
            Counter = Counter + 1

        end do

        Counter = 2

        do i = 1,Degree**2

            Node_Numberings(i,2) = Conversion_Vector(Counter)
            if (mod((real(Counter))/(Real(Degree)+1.0), 1.0) == 0) Counter = Counter + 1
            Counter = Counter + 1

        end do

        Counter = Degree + 3

        do i = 1,Degree**2

            Node_Numberings(i,3) = Conversion_Vector(Counter)
            if (mod((real(Counter)-(Real(Degree)+1.0))/(Real(Degree)+1.0), 1.0) == 0) Counter = Counter + 1
            Counter = Counter + 1

        end do

        Counter = Degree + 2

        do i = 1,Degree**2

            Node_Numberings(i,4) = Conversion_Vector(Counter)
            if (mod((real(Counter)-(Real(Degree)))/(Real(Degree)+1.0), 1.0) == 0) Counter = Counter + 1
            Counter = Counter + 1

        end do

    end subroutine Calculate_Node_Numberings

end module m_VTK_Writer