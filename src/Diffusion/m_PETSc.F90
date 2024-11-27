module m_PETSc
!
! Purpose:
! To solve the linear system Ax = b using PETSc
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ==========      =====================
! 31/10/2024    Daniel Beer     Original code
!
#include <petsc/finclude/petscksp.h>

use petscksp
use m_Read_Properties_D
implicit none

type :: PETScType
  private
  PetscScalar, pointer           :: xx_v(:)
  KSP                            :: ksp
  Mat                            :: A_CSR
  Vec                            :: x, b
  PC                             :: pc
end type

type :: PETScGroupType
  private
  integer, dimension(:), allocatable         :: indices
  type(PETScType), dimension(:), allocatable :: PETSc
end type

contains

  subroutine Start_PETSc(this,Properties,m,N,tol,max_iter,Periodic_Pairs)

    class(PETScGroupType)   :: this
    type(NTypeD), intent(in) :: N
    type(Mesh), intent(inout)  :: m
    type(PropertiesTypeD), intent(in)  :: Properties

    integer                                         :: ierr, Group_index, Nodes_per_Element, k, i, j
    double precision, intent(in)                    :: tol
    integer, intent(in)                             :: max_iter
    integer, dimension(:), allocatable              :: nnz
    integer, dimension(:,:), intent(in)             :: Periodic_Pairs
    real(kind = 8), dimension(2)                    :: Periodic_Value

    allocate(this%indices(N%N))
    allocate(this%PETSc(N%Group))
    allocate(nnz(N%N))

    call Calculate_NNZ(Properties, N, m, NNZ, Periodic_Pairs)

    do i = 1,N%N
      this%indices(i) = i-1
    end do

    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)

    Periodic_Value = [1e9,-1e9]

    do Group_index = 1,N%Group

      call VecCreate(PETSC_COMM_WORLD,this%PETSc(Group_index)%x,ierr)
      call VecSetSizes(this%PETSc(Group_index)%x, PETSC_DECIDE, N%N, ierr)
      call VecSetFromOptions(this%PETSc(Group_index)%x, ierr)
      call VecDuplicate(this%PETSc(Group_index)%x,this%PETSc(Group_index)%b,ierr)

      call MatCreateSeqAIJ(PETSC_COMM_SELF,N%N,N%N,PETSC_DECIDE,nnz,this%PETSc(Group_index)%A_CSR,ierr)
      call MatSetFromOptions(this%PETSc(Group_index)%A_CSR,ierr)
      call MatSetUp(this%PETSc(Group_index)%A_CSR,ierr)

      m%Cell_Pointers = m%Cell_Pointers - 1

      do k = 1, N%Element
        Nodes_per_Element = m%N_Cell_Pointers(k)
        do i = 1, Nodes_per_Element
            do j = 1, Nodes_per_Element
                call MatSetValues(this%PETSc(Group_index)%A_CSR, 1, m%Cell_Pointers(k,i), 1, m%Cell_Pointers(k,j), Properties%Elements(k)%K_Matrix(Group_index,i,j), ADD_VALUES, ierr)
            end do
        end do
      end do

      if (Properties%LBC == 4 .or. Properties%RBC == 4 .or. Properties%BBC == 4 .or. Properties%TBC == 4) then
        do i = 1, size(Periodic_Pairs,1)
          call MatSetValues(this%PETSc(Group_index)%A_CSR, 1, Periodic_Pairs(i,1), 1, Periodic_Pairs(i,1), Periodic_Value(1), ADD_VALUES, ierr)
          call MatSetValues(this%PETSc(Group_index)%A_CSR, 1, Periodic_Pairs(i,1), 1, Periodic_Pairs(i,2), Periodic_Value(2), ADD_VALUES, ierr)
          call MatSetValues(this%PETSc(Group_index)%A_CSR, 1, Periodic_Pairs(i,2), 1, Periodic_Pairs(i,1), Periodic_Value(2), ADD_VALUES, ierr)
          call MatSetValues(this%PETSc(Group_index)%A_CSR, 1, Periodic_Pairs(i,2), 1, Periodic_Pairs(i,2), Periodic_Value(1), ADD_VALUES, ierr)
        end do
      end if

      m%Cell_Pointers = m%Cell_Pointers + 1

      call MatAssemblyBegin(this%PETSc(Group_index)%A_CSR,MAT_FINAL_ASSEMBLY,ierr)
      call MatAssemblyEnd(this%PETSc(Group_index)%A_CSR,MAT_FINAL_ASSEMBLY,ierr)

      call KSPCreate(PETSC_COMM_WORLD,this%PETSc(Group_index)%ksp,ierr)
      call KSPSetOperators(this%PETSc(Group_index)%ksp,this%PETSc(Group_index)%A_CSR,this%PETSc(Group_index)%A_CSR,ierr)
      call KSPGetPC(this%PETSc(Group_index)%ksp,this%PETSc(Group_index)%pc,ierr)
      
      ! if (P_C == 0) then
      !   call PCSetType(this%PETSc(Group_index)%pc,PCNONE,ierr)
      ! else if (P_C == 1) then
      !   call PCSetType(this%PETSc(Group_index)%pc,PCJACOBI,ierr)
      ! else if (P_C == 2) then
      !   call PCSetType(this%PETSc(Group_index)%pc,PCCHOLESKY,ierr)
      ! else if (P_C == 3) then
      !   call PCSetType(this%PETSc(Group_index)%pc,PCLU,ierr)
      ! end if

      call KSPSetTolerances(this%PETSc(Group_index)%ksp,tol,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,max_iter,ierr)
      call KSPSetFromOptions(this%PETSc(Group_index)%ksp,ierr)

      ! call MatView(this%PETSc(Group_index)%A_CSR, PETSC_VIEWER_STDOUT_SELF, ierr)
      ! stop
 
    end do

  end subroutine 

  subroutine KSP_Solve(this,d,N,x_new,Group_index)

    class(PETScGroupType), intent(inout)          :: this
    double precision, dimension(:), intent(inout) :: x_new
    integer, intent(in)                           :: N
    double precision, dimension(:), intent(in)    :: d
    integer                                       :: ierr, its, Group_index

    call VecSetValues(this%PETSc(Group_index)%b, n, this%indices, d, INSERT_VALUES, ierr)
    call VecAssemblyBegin(this%PETSc(Group_index)%b, ierr)
    call VecAssemblyEnd(this%PETSc(Group_index)%b, ierr)
    ! call VecView(this%PETSc(Group_index)%b, PETSC_VIEWER_STDOUT_WORLD, ierr)

    call KSPSetType(this%PETSc(Group_index)%ksp, KSPCG, ierr)
    ! call KSPView(this%PETSc(Group_index)%ksp, PETSC_VIEWER_STDOUT_SELF, ierr)
    call KSPSolve(this%PETSc(Group_index)%ksp,this%PETSc(Group_index)%b,this%PETSc(Group_index)%x,ierr)
    call KSPGetIterationNumber(this%PETSc(Group_index)%ksp, its, ierr);

    !print *, 'Iterations: ', its

    call VecGetArrayF90(this%PETSc(Group_index)%x,this%PETSc(Group_index)%xx_v,ierr) 
    x_new = this%PETSc(Group_index)%xx_v
    call VecRestoreArrayF90(this%PETSc(Group_index)%x, this%PETSc(Group_index)%xx_v, ierr)

  end subroutine KSP_Solve

  subroutine End_PETSc(this)
    class(PETScGroupType), intent(inout) :: this

    integer :: ierr, i

    do i = 1,size(this%PETSc)

      call KSPDestroy(this%PETSc(i)%ksp,ierr)
      call MatDestroy(this%PETSc(i)%A_CSR,ierr)
      call VecDestroy(this%PETSc(i)%x,ierr)
      call VecDestroy(this%PETSc(i)%b,ierr)
      
    end do

    call PetscFinalize(ierr)

  end subroutine End_PETSc

  subroutine Calculate_NNZ(Properties, N, this, NNZ, Periodic_Pairs)

    type(PropertiesTypeD), intent(in) :: Properties
    type(NTypeD), intent(in)          :: N
    type(Mesh), intent(in)            :: this

    integer, dimension(:), intent(inout) :: NNZ

    integer, dimension(:,:), intent(in)  :: Periodic_Pairs

    logical, dimension(N%N) :: counter

    integer :: i, j, k, Nodes_per_Element

    do i = 1, N%N
        counter = .false.
        do j = 1, N%Element
            Nodes_per_Element = size(Properties%Elements(j)%Coordinates,1)
            if (any(this%Cell_Pointers(j,:) == i)) then
                do k = 1, Nodes_per_Element
                    counter(this%Cell_Pointers(j,k)) = .true.
                end do
            end if
        end do
        NNZ(i) = count(counter)
    end do

    if (Properties%LBC == 4 .or. Properties%RBC == 4 .or. Properties%BBC == 4 .or. Properties%TBC == 4) then
        do i = 1, size(Periodic_Pairs,1)
            NNZ(Periodic_Pairs(i,1)+1) = NNZ(Periodic_Pairs(i,1)+1) + 1
            NNZ(Periodic_Pairs(i,2)+1) = NNZ(Periodic_Pairs(i,2)+1) + 1
        end do
    end if

end subroutine Calculate_NNZ

end module m_PETSc