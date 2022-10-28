
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
  public :: factorial
  public :: unify_rows
  private :: sparse_gem_pig
  public :: block_jacobi_uncoupled
!  public :: sparse_gem_pig_AtA


    contains
    
  subroutine unify_rows(a,b)
       use typy
       use sparsematrix
       use global_objs
       use debug_tools

       !> matrix in sparse format
       class(smtx), intent(in out) :: a
       real(kind=rkind), dimension(:), intent(in out) :: b

       integer(kind=ikind) :: fin
       integer(kind=ikind), dimension(:), allocatable, save :: indices
       real(kind=rkind), dimension(:), allocatable, save :: values
       integer(kind=ikind) :: nelem, i, j
       real(kind=rkind) :: the_sum

       fin = a%getn()

       do i=1, fin
         call a%getrow(i=i, v=values, jj=indices, nelem=nelem)

         the_sum = sum(abs(values(1:nelem)))

         do j=1, nelem
           call a%set(values(j)/the_sum, i, indices(j))
         end do
         b(i) = b(i)/the_sum
       end do

     end subroutine unify_rows

    
    
        
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
      use core_tools


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
    



   subroutine block_jacobi(A,xvect,bvect,blindex, itcnt, repsexit, maxitcount, details)
      use typy
      use sparsematrix
      use pde_objs
      use solvers
      use debug_tools
      use globals
      use printtools
      use core_tools
      
      !> input global matrix
      class(smtx), intent(in out) :: A
      !> global vector with solution
      real(kind=rkind), dimension(:), intent(in out) :: xvect
      !> global b-vector
      real(kind=rkind), dimension(:), intent(in) :: bvect
      !> indeces of matrix diagonal blocks, 1st column start indices, 2nd column end indices
      integer(kind=ikind), dimension(:,:), intent(in) :: blindex
      !> iteration count 
      integer(kind=ikind), intent(out) :: itcnt
      !> required minimal residual
      real(kind=rkind) :: repsexit
      !> maximal iteration count
      integer(kind=ikind), intent(in) :: maxitcount
      !> if present and true prints Frobenius norms of submatrices
      logical, intent(in), optional :: details
      
      class(smtx), dimension(:,:), allocatable, save :: blockmat
      integer(kind=ikind) :: no_blocks
      integer(kind=ikind) :: finbig
      integer(kind=ikind), dimension(:), allocatable, save :: indices
      real(kind=rkind), dimension(:), allocatable, save :: values, xold, biter
      integer(kind=ikind) :: nelem, i, j, iblock, jblock, ilocal, jlocal, xlow, xhigh, blow, bhigh, pcg_it
      integer(kind=ikind) :: i_prev, j_prev, k
      real(kind=rkind) :: repstot, repslocal
      integer, save :: itfile
      real(kind=rkind), dimension(:,:), allocatable, save :: frobnorms
      character(len=4096) :: text
    
      
      finbig = A%getn()

      
      if (finbig /= A%getm()) then
        print *, "runtime error, matrix is not square matrix"
        print *, "exited from simplelinalg::block_jacobi"
        ERROR STOP
      end if
      
      
      if (ubound(xvect,1) /= finbig) then
        print *, "x vector has incorrect dimension"
        print *, "x vector has:", ubound(xvect,1), "components, matrix A has dimension:", finbig
        print *, "exited from simplelinalg::blockjacobi"
        ERROR STOP
      end if
      
      if (ubound(bvect,1) /= finbig) then
        print *, "b vector has incorrect dimension"
        print *, "b vector has:", ubound(bvect,1), "components, matrix A has dimension:", finbig
        print *, "exited from simplelinalg::blockjacobi"
        ERROR STOP
      end if      
      
      no_blocks = ubound(blindex,1)
      

      
      if (.not. allocated(blockmat)) then
        allocate(blockmat(no_blocks, no_blocks))
        do iblock=1, no_blocks
          do jblock=1, no_blocks
            call blockmat(iblock,jblock)%init(blindex(iblock,2)-blindex(iblock,1)+1, blindex(iblock,2)-blindex(iblock,1)+1)
          end do
        end do
        allocate(xold(ubound(xvect,1)))
        allocate(biter(ubound(xvect,1)))
        open(newunit=itfile, file="out/BJ.iterations", status="replace", action="write")
        call print_logo(itfile)
        if (.not. present(details)) then
          write(unit=itfile, fmt=*) "# time                            iterations       final iteration increment" 
        else
          write(unit=itfile, fmt=*) "# time                            iterations       final iteration increment,   frobnorm(1,2) &
                    .          frobnorm(2,1)"
          write(unit=itfile, fmt=*) "#                                                                               ------------- &
                    .          -------------"
          write(unit=itfile, fmt=*) "#                                                                               frobnorm(1,1) &
                    .         frobnorm(2,2)"
          write(unit=itfile, fmt=*) "# "          
          write(unit=itfile, fmt=*) "#  ------------------------------------------------------------------------------------------"
        end if
      end if
      
      
      do i=1, finbig
        call a%getrow(i=i, v=values, jj=indices, nelem=nelem)
        do j=1, nelem
          do k=1, no_blocks
            if (i >= blindex(k,1) .and. i <= blindex(k,2) ) then
              iblock = k
            end if
            
            if (indices(j) >= blindex(k,1) .and. indices(j) <= blindex(k,2) ) then
              jblock = k
            end if
          end do
        
          
          if (iblock > 1) then
            i_prev = blindex(iblock-1,2)
          else
            i_prev = 0
          end if
          
          if (jblock > 1) then
            j_prev = blindex(jblock-1,2)
          else
            j_prev = 0
          end if
            
          
          ilocal = i - i_prev
          jlocal = indices(j) - j_prev
          

          call blockmat(iblock,jblock)%set(values(j), ilocal, jlocal)


        end do
      end do
          

      itcnt = 0
      
      if (present(details)) then
        if (details) then
          if (.not. allocated(frobnorms)) allocate(frobnorms(ubound(blockmat,1),ubound(blockmat,2) ))
          frobnorms = 0
          do iblock=1, ubound(blockmat,1)
            do jblock=1, ubound(blockmat,2)
              do i=1, blockmat(iblock, jblock)%getn()
                call blockmat(iblock, jblock)%getrow(i=i, v=values, jj=indices, nelem=nelem)
                do j=1, nelem
                  frobnorms(iblock, jblock) = frobnorms(iblock, jblock) + abs(values(j))
                end do
              end do
            end do  
          end do        
        end if
      end if

      do iblock=1,  ubound(blockmat,1)
        call LDUd(blockmat(iblock, iblock))
      end do

	
      do    
        itcnt = itcnt + 1
        xold = xvect
        do iblock=1, ubound(blockmat,1)
          blow = blindex(iblock,1)
          bhigh = blindex(iblock,2)
          do jblock=1, ubound(blockmat,1)
            xlow = blindex(jblock,1)
            xhigh = blindex(jblock,2)
            if (iblock /= jblock) then          
              biter(blow:bhigh) = bvect(blow:bhigh) - blockmat(iblock, jblock)%mul(xold(xlow:xhigh))
            end if
          end do
        end do
        
        do iblock=1, ubound(blockmat,1)
          blow = blindex(iblock,1)
          bhigh = blindex(iblock,2)
          call LDUback(blockmat(iblock, iblock), biter(blow:bhigh), xvect(blow:bhigh))
        end do  
  
        repstot = norm2(xvect - xold)
       
        if (repstot < repsexit) then
          write(unit=terminal, fmt=*) "Block Jacobi iteration count:", itcnt
          write(unit=itfile, fmt=*) time, itcnt, repstot, frobnorms(1,2)/frobnorms(1,1), frobnorms(2,1)/frobnorms(2,2)
          EXIT
        end if
        
        if (itcnt > maxitcount) then
          call write_log(text="Block Jacobi failed to converge, the maximal allowed number of iterations was:", int1=maxitcount)
          if (present(details)) then
            if (details) then
              write(unit=text, fmt=*) frobnorms(1,2)/frobnorms(1,1), frobnorms(2,1)/frobnorms(2,2)
              call write_log(text="Frobenius submatrices norms were:", text2=cut(text))
            end if
          end if
              
          ERROR STOP
        end if
        
      end do
      
	end subroutine block_jacobi
                  
    
    subroutine block_jacobi_uncoupled(A,xvect,bvect,blindex, itcnt, repsexit, maxitcount)
      use typy
      use sparsematrix
      use pde_objs
      use solvers
      use debug_tools
      use globals
      use global4solver
      use printtools
      use core_tools
      
      !> input global matrix
      class(smtx), intent(in out) :: A
      !> global vector with solution
      real(kind=rkind), dimension(:), intent(in out) :: xvect
      !> global b-vector
      real(kind=rkind), dimension(:), intent(in) :: bvect
      !> indeces of matrix diagonal blocks, 1st column start indices, 2nd column end indices
      integer(kind=ikind), dimension(:,:), intent(in) :: blindex
      !> iteration count for PCG applied on submatrix (2,2)
      integer(kind=ikind), intent(out) :: itcnt
      !> required minimal residual for PCG applied on submatrix (2,2)
      real(kind=rkind) :: repsexit
      !> maximal iteration count for PCG applied on submatrix (2,2)
      integer(kind=ikind), intent(in) :: maxitcount

     
      
      class(extsmtx), dimension(:,:), allocatable, save :: blockmat
      integer(kind=ikind) :: no_blocks
      integer(kind=ikind) :: finbig
      integer(kind=ikind), dimension(:), allocatable, save :: indices
      real(kind=rkind), dimension(:), allocatable, save :: values
      integer(kind=ikind) :: nelem, i, j, iblock, jblock, ilocal, jlocal, xlow, xhigh, blow, bhigh, pcg_it
      integer(kind=ikind) :: i_prev, j_prev, k
      real(kind=rkind) :: repstot, repslocal
      integer, save :: itfile
      real(kind=rkind), dimension(:,:), allocatable, save :: frobnorms
      character(len=4096) :: text
      
      integer(kind=ikind), dimension(:), allocatable, save :: p1, p2
      logical, save :: solved = .false.
      real(kind=rkind), dimension(:), allocatable :: blocal
    
      finbig = A%getn()

      
      if (finbig /= A%getm()) then
        print *, "runtime error, matrix is not square matrix"
        print *, "exited from simplelinalg::block_jacobi"
        ERROR STOP
      end if
      
      
      if (ubound(xvect,1) /= finbig) then
        print *, "x vector has incorrect dimension"
        print *, "x vector has:", ubound(xvect,1), "components, matrix A has dimension:", finbig
        print *, "exited from simplelinalg::blockjacobi"
        ERROR STOP
      end if
      
      if (ubound(bvect,1) /= finbig) then
        print *, "b vector has incorrect dimension"
        print *, "b vector has:", ubound(bvect,1), "components, matrix A has dimension:", finbig
        print *, "exited from simplelinalg::blockjacobi"
        ERROR STOP
      end if      
      
      no_blocks = ubound(blindex,1)
      
      
      if (.not. allocated(blockmat)) then
        allocate(blockmat(no_blocks, no_blocks))
        do iblock=1, no_blocks
          do jblock=1, no_blocks
            call blockmat(iblock,jblock)%init(blindex(iblock,2)-blindex(iblock,1)+1, blindex(iblock,2)-blindex(iblock,1)+1)
          end do
        end do
      end if
      
      
      do i=1, finbig
        call a%getrow(i=i, v=values, jj=indices, nelem=nelem)
        do j=1, nelem
          do k=1, no_blocks
            if (i >= blindex(k,1) .and. i <= blindex(k,2) ) then
              iblock = k
            end if
            
            if (indices(j) >= blindex(k,1) .and. indices(j) <= blindex(k,2) ) then
              jblock = k
            end if
          end do
        
          
          if (iblock > 1) then
            i_prev = blindex(iblock-1,2)
          else
            i_prev = 0
          end if
          
          if (jblock > 1) then
            j_prev = blindex(jblock-1,2)
          else
            j_prev = 0
          end if
            
          
          ilocal = i - i_prev
          jlocal = indices(j) - j_prev
          

          call blockmat(iblock,jblock)%set(values(j), ilocal, jlocal)


        end do
      end do
          
      

	  if (.not. allocated(p1)) then
		allocate(p1(blindex(1,2)))
		allocate(p2(blindex(1,2)))
	  end if

	  call LDUd(blockmat(1,1), pivtype=5, ilev=0, perm1=p1, perm2=p2)
	  call LDUback(blockmat(1,1), bvect(blindex(1,1):blindex(1,2)), xvect(blindex(1,1):blindex(1,2)), p1=p1, p2=p2)
	  
	  if (.not. solved) itcnt = 0

	 
	  if (solved) then 
		if (solver_name == "LDUPCG") then
          call diag_precond(a=blockmat(2,2), x=xvect(blindex(2,1):blindex(2,2)), mode=1)
          call CGnormal(A=blockmat(2,2), b=bvect(blindex(2,1):blindex(2,2)),x=xvect(blindex(2,1):blindex(2,2)), &
						ilev1=0,itmax1=maxitcount,reps1=repsexit, itfin1=itcnt)
          call diag_precond(a=blockmat(2,2), x=xvect(blindex(2,1):blindex(2,2)), mode=-1)
         else if (solver_name == "LDULDU") then
			call LDUd(blockmat(2,2),  pivtype=5, perm1=p1, perm2=p2)
			call LDUback(blockmat(2,2), bvect(blindex(2,1):blindex(2,2)), xvect(blindex(2,1):blindex(2,2)), p1=p1, p2=p2) 
			itcnt = 0
		 else
			STOP "exited from simplelingalg::block_jacobi_uncoupled"
		 end if
	 end if
	 
	 solved = .true.

    
    end subroutine block_jacobi_uncoupled
    
    
    
   
    
    
      !> this procedure bears a codename sparse_gem_pig because it creates full matrix out of sparse matrix and solves it on using full Gauss elimination, only for debugging
   subroutine sparse_gem_pig(A,b,x)!,ptype,ilev,ierr, itmax, reps_rel, reps_abs, info, it)
     use linalg
     use typy
     use sparsematrix
     
     class(smtx), intent(in out) :: a
     real(kind=rkind), dimension(:), intent(in) :: b
     real(kind=rkind), dimension(:), intent(in out) :: x
!   !> level of information, default = 0
!       integer, intent(in), optional                  :: ilev
!       !> kind of preconditioner, default 0
!       integer, intent(in), optional                  :: ptype
!       !> error message\n 0 .. OK\n
!       !>                1 .. after itmax iterions is not founded sufficiently small relative error
!       integer, intent(out), optional                 :: ierr
!       !> maximum allowed iterations, default = 500
!       integer(kind=ikind), intent(in), optional      :: itmax
!       !> maximal relative error
!       real(kind=rkind), intent(in), optional         :: reps_rel
!       !> maximal absolute error
!       real(kind=rkind), intent(in), optional         :: reps_abs
!       !> operations count
!       type(info_type), intent(inout), optional       :: info
!       !> iteration number
!       integer(kind=ikind), intent(out), optional     :: it
     real(kind=rkind), dimension(:,:), allocatable :: matice
     integer(kind=ikind), dimension(:), allocatable :: jj
     real(kind=rkind), dimension(:), allocatable :: v
     integer(kind=ikind) :: nelem
     integer(kind=ikind) :: i, j
     
     allocate(matice(ubound(b,1), ubound(b,1)))
 
     matice = 0.0_rkind
 
     do i=1,a%getn()
       call a%getrow(i,v,jj,nelem)
       do j=1, nelem
       	 matice(i, jj(j)) = v(j)
       end do	
     end do
 
     call gem(matice, b, x)
 
   end subroutine sparse_gem_pig
   
   
   subroutine Ax4gmres(b,x,nsize)
	 use typy
	 use globals
	 
	 integer(kind=ikind), intent (in):: nsize
	 real(kind=rkind), dimension(1:nsize), intent(in) :: x
	 real(kind=rkind), dimension(1:nsize), intent(inout) :: b
	 
	 b = spmatrix%mul(x)
	 
   end subroutine Ax4gmres
         
         
   subroutine dummycond4gmres(b,x,nsize)
	use typy
    integer(kind=ikind), intent (in):: nsize
    real(kind=rkind), dimension(1:nsize), intent(in) :: x
    real(kind=rkind), dimension(1:nsize), intent(inout) :: b
    
    b=x
    
   end subroutine dummycond4gmres       
! 
!   
!     !> this procedure bears a codename sparse_gem_pig because it creates full matrix out of sparse matrix and solves it on using full Gauss elimination, only for debugging
!     subroutine sparse_gem_pig_AtA(A,b,x,ptype,ilev,ierr, itmax, reps, info,normmul)
!         use mytypes
!         use typy
!         use linAlg
!         implicit none
!         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         ! parametry
!         !> system matrix - supposed be symmetric positive definite\n
!         !> inout because of possibility of initalizing of preconditioner
!         type(smtx), intent(in out)                  :: A
!         !> right hand side
!         real(kind=rkind), dimension(:), intent(in) :: b
!         !> solution, on input initial aproximation
!         real(kind=rkind), dimension(:), intent(in out) :: x
!         !> level of information, default = 0
!         integer, intent(in), optional           :: ilev
!         !> kind of preconditioner, default 0
!         integer, intent(in), optional           :: ptype
!         !> error message\n 0 .. OK\n
!         !>                1 .. after itmax iterions is not founded sufficiently small relative error
!         integer, intent(out), optional          :: ierr
!         !> maximum allowed iterations, default = 500
!         integer(kind=ikind), intent(in), optional           :: itmax
!         !> wanted relative error, default=  1e-5
!         real(kind=rkind), intent(in), optional  :: reps
!         !> operations count
!         type(info_type), intent(inout), optional  :: info
!         !> zda pro normeq=true delat prenasobeni
!         logical, optional, intent(in) :: normmul
!     real(kind=rkind), dimension(:,:), allocatable :: matice, matice2
!     real(kind=rkind), dimension(:), allocatable :: Atb
!     integer(kind=ikind) :: i
!     
!     allocate(matice(ubound(b,1), ubound(x,1)))
!     allocate(matice2(ubound(x,1), ubound(b,1)))
!     allocate(Atb(ubound(x,1)))
!     
!     matice = 0.0_rkind
!     matice2 = 0.0_rkind
! 
!     do i=1,ubound(a%vals,1)
!       if (a%ii(i) > 0 .or. a%jj(i) > 0) then
!         matice(a%ii(i), a%jj(i)) =  matice(a%ii(i), a%jj(i)) + a%vals(i)
!       end if
!     end do
! 
!     do i=1,ubound(a%vals,1)
!       if (a%ii(i) > 0 .or. a%jj(i) > 0 .and.  abs(a%vals(i)) > epsilon(a%vals(i))) then
!         matice2(a%jj(i), a%ii(i)) =  matice2(a%jj(i), a%ii(i)) + a%vals(i)
!       end if
!     end do
!     
! 
!     
!     matice = matmul(matice2,matice)
! 
!     
!     Atb = matmul(matice2,b)
!     
! !     call printmtx(matice)
!     
! !     call printmtx(Atb)
!     
!     call gem(matice, Atb, x)
!   end subroutine sparse_gem_pig_AtA
 

end module simplelinalg
