
! Copyright 2008 Michal Kuraz, Petr Mayer, Copyright 2016  Michal Kuraz, Petr Mayer, Johanna Bloecher


! This file is part of DRUtES.
! DRUtES is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! DRUtES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with DRUtES. If not, see <http://www.gnu.org/licenses/>.


!> \file simplelinalg.f90
!! \brief Diagonal preconditioner -- simple but robust for diagonaly dominant problems.
!<





module simplelinalg
  public :: diag_precond
  public :: invert_matrix
  private ::  determinant
  public :: factorial


    contains
        
          !> right hand side diagonal preconditioner - preprocessor and postprocesor
    !! preforms the diagonal matrix scaling as described in Kuraz & Mayer: Algorithms for solving Darcian flow in structured porous media
    !! 
    !! the postprocesor performs following operation with the vector of solution
    !! \f[ \mathbf{x} = \mathbf{x} \times \frac{1}{\mathbf{a_{diag}}} \f]
    !! where a_diag is a vector of diagonal values in solution matrix \n
    !! parameter mode defines the preprocessing or postprocessing mode for this routine   
    !>
    subroutine diag_precond(a, x, prmt, mode)
      use typy
      use sparsematrix
      use globals
      use globals2D
      use debug_tools

      !> matrix in sparse format
      class(extsmtx), intent(in out) :: a
      !> the solution vector to be passed into zero iteration, it is formed out of the solution in previous time step level
      real(kind=rkind), dimension(:), intent(in out), optional :: x
      !> permut vector, usefull for domain decomposition
      integer(kind=ikind), dimension(:), intent(in), optional :: prmt
      !> mode defines whether to procedure should be stared as a preprocesor or postprocesor \n
      !! mode = 1  -> preprocessor
      !! mode = -1 -> postprocesor
      !<
      integer, intent(in) :: mode
      integer(kind=ikind), dimension(:), allocatable, save :: indexes
      integer(kind=ikind) :: i, finn, finm, nelem, j, k, diag_row
      real(kind=rkind), dimension(:), allocatable, save :: values
      real(kind=rkind) :: tmp
      
      finn=a%getn()
      finm=a%getm()

      if (.not. allocated(indexes)) then
        allocate(indexes(2*drutes_config%dimen+1))
        allocate(values(2*drutes_config%dimen+1))
      end if


      if (.not. allocated(a%weight)) then
        allocate(a%weight(finm))
      end if
      
      select case(mode)
        case(1)
            a%weighted = .true.
            do i=1, finm
              if (present(prmt)) then
                diag_row = prmt(i)
              else
                diag_row = i
              end if
              a%weight(i) = 1.0_rkind/a%get(diag_row,i)
            end do

            do i=1, finn
              call a%getrow(i=i, v=values, jj=indexes, nelem=nelem)
              do k=1, nelem
                j = indexes(k)
                tmp = a%get(i,j)
                call a%set(tmp*a%weight(j),i,j)
              end do
            end do
    
            if (present(x)) then
              x(1:finm) = x(1:finm)/a%weight(1:finm)
            end if
  

        case(-1)
          if (.not. a%weighted) then
            error stop "ERROR the matrix was not priorly diagonalized, called from simplelinalg::diag_precond"
          end if

          if (.not. present(x)) then
            print *, "x-vector not present, incorrect function call, called from simplelinalg::diag_precond"
            error stop
          end if

          x(1:finm)=x(1:finm)*a%weight

          do i=1, finn
            call a%getrow(i=i, v=values, jj=indexes, nelem=nelem)
            do k=1, nelem
              j = indexes(k)
              tmp = a%get(i,j)
              call a%set(tmp/a%weight(j),i,j)
            end do
          end do
          a%weighted = .false.

        case default
          print *, "RUNTIME error in parameter mode passed into procedure fem_tools::diag_precond"
          ERROR STOP
      end select


    end subroutine diag_precond

    !> subroutine, that evaluates inverse matrix up to dimension 3
    subroutine invert_matrix(A)
      use typy


      real(kind=rkind), dimension(:,:), intent(in out) :: A
      real(kind=rkind), dimension(3,3) :: A_loc
      integer(kind=ikind) :: i,j, n

      n = ubound(A,1)

      select case(n)
        case(1)
            A(1,1) = 1.0_rkind/A(1,1)
        case(2)
            A_loc(1:n,1:n) = A
            A(1,1) = 1.0_rkind/determinant(A_loc(1:2,1:2))*A_loc(2,2)
            A(1,2) = -1.0_rkind/determinant(A_loc(1:2,1:2))*A_loc(1,2)
            A(2,1) = -1.0_rkind/determinant(A_loc(1:2,1:2))*A_loc(2,1)
            A(2,2) =  1.0_rkind/determinant(A_loc(1:2,1:2))*A_loc(1,1)
        case(3)
            print *, "no implemented, terminated from fem_tools::invert_matrix"
            ERROR stop
      end select

    end subroutine invert_matrix

    !> the matrix dimension must be exactly (2,2) or (3,3)
    function determinant(A) result(det)
      use typy


      real(kind=rkind) :: det
      real(kind=rkind), dimension(:,:), intent(in) :: A


      if (ubound(A,1) /= ubound(A,2)) then
        print *, "ERROR: input matrix is not a square matrix, called from fem_tools::determinant()"
        ERROR STOP
      end if


      select case(ubound(A,1))
        case(2)
              det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
        case(3)
              det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2)&
                  - A(1,3)*A(2,2)*A(3,1) - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)
        case default
              print *, "ERROR: incorrect matrix dimension, called from fem_tools::determinant()"
              ERROR STOP
      end select


    end function determinant
    
    function factorial(in) result(fact)
      use typy
      
      integer(kind=ikind), intent(in) :: in
      integer(kind=ikind) :: fact
      
      integer(kind=ikind) :: i
      
      fact = in
      do i=in-1, 1, -1
        fact = fact*i
      end do
      
    end function factorial

 

end module simplelinalg
