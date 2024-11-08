module m_Small_Matrix_Solver
!
! Purpose:
! To solve small matrices using the analytical solution or LAPACK.
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
  implicit none

  contains

    subroutine Solve_Matrix(A, b, x)
      real(8), dimension(:, :), intent(in) :: A
      real(8), dimension(size(A,1), size(A,2)) :: A_new
      real(8), dimension(:), intent(in) :: b
      real(8), dimension(:), intent(out) :: x
      real(8) :: det
      integer, dimension(size(A,1)) :: IPIV
      integer :: INFO

      if (size(A, 1) == 2 .and. size(A, 2) == 2) then    

        det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)

        x(1) = (A(2, 2) * b(1) - A(1, 2) * b(2)) / det
        x(2) = (-A(2, 1) * b(1) + A(1, 1) * b(2)) / det

      else if (size(A, 1) == 3 .and. size(A, 2) == 3) then

        det = A(1, 1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
            - A(1, 2) * (A(2, 1) * A(3, 3) - A(2, 3) * A(3, 1)) &
            + A(1, 3) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1))

        x(1) = (b(1) * (A(2, 2) * A(3, 3) - A(2, 3) * A(3, 2)) &
             - b(2) * (A(1, 2) * A(3, 3) - A(1, 3) * A(3, 2)) &
             + b(3) * (A(1, 2) * A(2, 3) - A(1, 3) * A(2, 2))) &
             / det

        x(2) = (b(1) * (A(2, 3) * A(3, 1) - A(2, 1) * A(3, 3)) &
             - b(2) * (A(1, 3) * A(3, 1) - A(1, 1) * A(3, 3)) &
             + b(3) * (A(1, 3) * A(2, 1) - A(1, 1) * A(2, 3))) &
             / det

        x(3) = (b(1) * (A(2, 1) * A(3, 2) - A(2, 2) * A(3, 1)) &
              - b(2) * (A(1, 1) * A(3, 2) - A(1, 2) * A(3, 1)) &
              + b(3) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1))) &
              / det

      else

        x = b

        A_new = A

        call dgesv(size(A, 1), 1, A_new, size(A, 1), IPIV, x, size(A, 1), INFO)

      end if
      
    end subroutine Solve_Matrix

end module m_Small_Matrix_Solver
