module femmat
  public :: solve_picard
  private :: results_extractor
  public :: assemble_mat
  

  contains
  !>solves the problem using Picard method
    subroutine solve_picard(ierr, itcount, success)
      use typy
      use globals
      use global_objs
      use pde_objs
      use RE_constitutive
      use feminittools
      use fem_tools
      use linAlg
      use debug_tools
      use solvers
      use postpro
      use decomp_tools
      use decomp_vars
      use simplelinalg
      use debug_tools
      use mtxiotools
!       use decomposer


      integer, intent(out) :: ierr
      integer(kind=ikind), intent(out) :: itcount
      integer(kind=ikind) ::  fin, pcg_it
      logical, intent(out) :: success
      integer(kind=ikind) :: i, proc,j, m, l, k, top, bot
      real(kind=rkind) :: error, loc_error
      logical :: dt_fine
      integer :: ierr_loc
      real(kind=rkind), dimension(:), allocatable :: vcttmp
      real(kind=rkind) :: lambda_l, lambda_h

      
      proc = ubound(pde,1)
      fin = maxval(pde(proc)%permut(:))
      itcount = 0
      allocate(vcttmp(ubound(pde_common%bvect,1)))
      

      do

	do proc=1, ubound(pde,1)
	  call icond4neumann(proc)
	end do

	itcount = itcount + 1
	
	call assemble_mat(ierr)
	
	
	
! 	call estimeigvalues(spmatrix, lambda_l, lambda_h)
! 	
! 	print *, lambda_l, lambda_h
	
	if (drutes_config%dimen >  0) then
	  call diag_precond(a=spmatrix, x=pde_common%xvect(1:fin,3), mode=1)
	end if

	call solve_matrix(spmatrix, pde_common%bvect(1:fin), pde_common%xvect(1:fin,3),  itmax1=fin, &
		  reps1=1e-15_rkind)
! 

	if (drutes_config%dimen >  0) then
	  write(unit=file_itcg, fmt = *) time, pcg_it, itcount
	  call flush(file_itcg)
	  call diag_precond(a=spmatrix, x=pde_common%xvect(1:fin,3), mode=-1)
	end if
	

         error = norm2(pde_common%xvect(1:fin,2)-pde_common%xvect(1:fin,3))/ubound(pde_common%xvect,1)

         
	if (itcount == 1 .or. error <= iter_criterion) then
	  do proc=1, ubound(pde,1)
	    top = pde(proc)%permut(1)
	    bot = pde(proc)%permut(ubound(pde(proc)%permut,1))
	    dt_fine = pde(proc)%dt_check()
	  end do
        else
          dt_fine = .true.
        end if


	if (.not.(dt_fine)) then
	  ierr = 2
	  success = .false.
	  RETURN
	end if



	if (error <= iter_criterion) then
	  
	
	  ierr = 0
	  call results_extractor()
	  

	  success = .true.
	  RETURN
	else
	  pde_common%xvect(:,2) = pde_common%xvect(:,3)
	end if

	if (itcount >= max_itcount) then
	  ierr = 1
	  success = .false.
	  EXIT
	end if
      end do

    end subroutine solve_picard
    
    subroutine assemble_mat(ierr)
      use typy
      use globals
      use global_objs
      use pde_objs
      use capmat
      use stiffmat
      use feminittools
      use geom_tools
      use fem_tools
      use re_constitutive
      use linAlg
      use solver_interfaces
      use debug_tools

      integer, intent(out) :: ierr

      integer(kind=ikind) :: i,j,k,l, proc, ll, limits
            
      type(integpnt_str) :: quadpnt
      
      real(kind=rkind) :: value

      call null_problem(spmatrix, pde_common%bvect)
      
      limits = ubound(stiff_mat,1)/ubound(pde,1)



      do i=1, elements%kolik
	pde_common%current_el = i
	processes: do proc=1, ubound(pde,1)
		    do j=1+(proc-1)*limits, ubound(elements%data,2) + (proc-1)*limits
		      ll = j - (proc-1)*ubound(stiff_mat,1)/ubound(pde,1)
		      k = pde(proc)%permut(elements%data(i,ll))
		      if (k > 0) then
			elnode_prev(j) = pde_common%xvect(k,1)
		      else
			k = nodes%edge(elements%data(i,ll))
			call pde(proc)%bc(k)%value_fnc(pde(proc), i, ll, value)
			elnode_prev(j) = value
		      end if
		    end do
	end do processes
	

	call build_bvect(i, time_step)
	
	call build_stiff_np(i, time_step)
	
	call pde_common%time_integ(i)
		
	stiff_mat = stiff_mat + cap_mat


	call in2global(i,spmatrix, pde_common%bvect)

      end do


    end subroutine assemble_mat
    
    
    subroutine results_extractor()
      use typy
      use globals
      use global_objs
      use pde_objs

      integer(kind=ikind) :: i, j, proc, nd, edge
      real(kind=rkind) :: value

      
      do proc=1, ubound(pde,1)
	do i=1, elements%kolik
	  do j=1, ubound(elements%data,2)
	    nd = elements%data(i,j)
	    if (pde(proc)%permut(nd) > 0) then
	      pde(proc)%solution(nd) = pde_common%xvect(pde(proc)%permut(nd),3)
	    else
	      edge = nodes%edge(nd)
	      call pde(proc)%bc(edge)%value_fnc(pde(proc), i, j, value)
	      pde(proc)%solution(nd) = value
	    end if
	  end do
	end do
      end do

      pde_common%xvect(:,1) = pde_common%xvect(:,3)
      pde_common%xvect(:,2) = pde_common%xvect(:,3)

    end subroutine results_extractor
    


end module femmat