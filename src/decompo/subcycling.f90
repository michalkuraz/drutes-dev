module subcycling
    public :: subcyc_schwarz
    private :: build_xvect
    private :: domains_solved
    private :: results_extractor
    private :: update_domain
    private :: correctdoms
    private :: getlocres
! 
    contains

      !>solves the problem using Picard method with Schwarz multiplicative domain decomposition reflecting the nonlinear problem
      subroutine subcyc_schwarz(ierr, itcount, success)
	use typy
	use globals
	use global_objs
	use pde_objs
	use feminittools
	use fem_tools
	use linAlg
	use debug_tools
	use sparsematrix
	use decomp_vars
	use decomp_tools
	use decomposer
	use postpro
        use femmat
        use core_tools
        use simplelinalg
        use printtools
        use debug_tools
        use schwarz_dd

	integer, intent(out) :: ierr
	integer(kind=ikind), intent(out) :: itcount
	integer(kind=ikind) ::  fin, pcg_it, subfin
	logical, intent(out) :: success
	integer(kind=ikind) :: i, proc,j, m, l, k, ii, sub
	real(kind=rkind) :: error, loc_error, res_error, reps
	logical :: dt_fine, first_run, picard_done, reset_domains, bigupdate
	integer :: ierr_loc
	real(kind=rkind), dimension(:), allocatable, save :: resvct, corrvct
	character(len=128) :: text
	real(4), dimension(2,2) :: taaray
	real(4), dimension(2) :: rslt

! 	

        call etime(taaray(1,:), rslt(1))
	proc = ubound(pde,1)
	fin = maxval(pde(proc)%permut(:))
	itcount = 0
	inner_criterion = 1e-20

	!reset local-local cluster iteration count
        ddcoarse_mesh(:)%iter_count = 1
	
	!the residual vector is not allocated, thus we are at the beginning
	if (.not. allocated(resvct)) then
          first_run = .true.
	  allocate(resvct(fin))
	  allocate(corrvct(fin))          
          call set_subdomains()
          subdomain(:)%time = 0
          else
          first_run = .false.
        end if

        do i=1, ubound(subdomain,1)
	  if (subdomain(i)%critical) then 
	    subdomain(i)%time_step = time_step/3
	  else
	    subdomain(i)%time_step = time_step
	  end if
	end do

                
        !reset local subdomain iteration count
        subdomain(:)%itcount = 0                

        !reset successes on domains
        subdomain(:)%solved = .false.
   	
	picard_done = .false.
	

        write(unit=terminal, fmt="(a)") "  "
        write(unit=terminal, fmt="(a, I4,a)") " Solving", ubound(subdomain,1),  " subdomains .... "
        
        !compute the local residuum
	!and create the local matrices
	do i=1, ubound(subdomain,1)
	  subdomain(i)%xvect(:,1:3) = pde_common%xvect(subdomain(i)%permut(1:subdomain(i)%ndof),1:3)
	  if (subdomain(i)%critical) then
	    subdomain(i)%time_step = time_step/3
	  end if
	end do

        
        time_loop: do
	  
	  do sub=1, ubound(subdomain,1)
	    
	    call update_domain(sub, bigupdate)
	    
	    if (bigupdate) then
	      call correctdoms()
	    end if
	    
	    if (subdomain(sub)%time < time + time_step) then
	      
	      itcount = 0
	      
	      !solve nonlinearities on subdomain
	      picard: do
	      
		call locmat_assembler(subdomain(i), ierr, i)
	      
		itcount = itcount + 1
		
		call getlocres(subdomain(sub), resvct)
	        
	        corrvct(1:subfin) = 0.0

		if (subdomain(i)%critical) then
		  reps = 1e-20
		else
		  reps = 1e-10
		end if
		
		call diag_precond(a=subdomain(i)%matrix, prmt=subdomain(i)%permut,  mode=1)
		
		call solve_matrix(subdomain(i)%matrix, resvct, corrvct(1:subfin), ilev1=0, itmax1=subfin, reps1=reps)

		call diag_precond(a=subdomain(i)%matrix, x=corrvct(1:subfin), mode=-1)	
		
		error = maxval(abs(corrvct(1:subfin)))
		
		subdomain(i)%xvect(:,3) = subdomain(i)%xvect(:,2) + corrvct(1:subfin)
		
		call search_error_cluster(subdomain(i), itcount) 
		
		subdomain(i)%itcount = itcount
		
		if (error <= iter_criterion) then
		  subdomain(sub)%time = subdomain(sub)%time + time_step
		  EXIT picard
		else if (itcount > max_itcount) then
		  subdomain(sub)%time_step = subdomain(sub)%time_step*0.5
		end if
	      end do picard
	    end if
	  end do
	      
	  if (domains_solved() == ubound(subdomain,1)) then
	    call etime(taaray(2,:), rslt(2))
	    EXIT time_loop
	 end if
	end do time_loop
	    
		  
	call build_xvect()

	    
	do i=1, ubound(subdomain,1)
	  if (.not. subdomain(i)%critical .and. subdomain(i)%itcount>1 ) then
	    reset_domains = .true.
	    EXIT
	  else
	    reset_domains = .false.
	  end if
	end do
	      
        call printtime("Time spent on iterations: ", rslt(2)-rslt(1))
	if (reset_domains) call set_subdomains()
	success = .TRUE.
	ierr = 0
	write(unit=file_ddinfo, fmt=*) time, itcount, "|", subdomain(:)%ndof
	call flush(file_ddinfo)
	     


      end subroutine subcyc_schwarz
      
      subroutine build_xvect()
	use typy
	use globals
	use global_objs
	use pde_objs
	use decomp_vars
	use debug_tools
	
	integer(kind=ikind) :: i, glob, loc, subcrit


        pde_common%xvect(:,2:3) = 0.0_rkind
        
        
        do i = 1, ubound(subdomain,1)
	  do loc = 1, subdomain(i)%ndof
	    glob =  subdomain(i)%permut(loc)
	    pde_common%xvect(glob,2:3) = pde_common%xvect(glob,2:3) + subdomain(i)%xvect(loc,2:3)*prolong_mtx%get(glob,i)
	  end do
	end do
	

      end subroutine build_xvect
      
      subroutine update_domain(id, bigupdate)
	use typy
	use decomp_vars
	
	integer(kind=ikind), intent(in) :: id
	logical, intent(out) :: bigupdate
	integer(kind=ikind), dimension(:), allocatable, save :: nzindexes
	integer(kind=ikind) :: iglob, num, i, sub, j, iloc
	real(kind=rkind), dimension(:), allocatable, save :: weight
	real(kind=rkind) :: value, subtime, currval
	
	if (.not. allocated(nzindexes)) then
	  allocate(nzindexes(ubound(subdomain,1)))
	  allocate(weight(ubound(subdomain,1)))
	end if
	
	subtime = subdomain(id)%time
	bigupdate = .false.
	do i=1, ubound(subdomain(id)%xvect,1)
	  iglob = subdomain(id)%permut(i)
	  call prolong_mtx%getrow(i=iglob, v=weight, jj=nzindexes, nelem=num)
	  if (num > 1) then
	    currval = subdomain(id)%xvect(i,2)
	    value = 0
	    do j=1, num
	      sub = nzindexes(j)
	      iloc = subdomain(sub)%invpermut(iglob)
	      value = value + subdomain(sub)%returnval(iloc, subtime)*weight(j)
	    end do
	    if (abs(currval - value) > iter_criterion) then
	      bigupdate = .true.
	    end if
	    subdomain(id)%xvect(i, 2:3) = value
	  end if
	end do
	  
      
      
      end subroutine update_domain
      


      subroutine search_error_cluster(sub, itcount)
	use typy
	use globals
	use global_objs
	use decomp_vars
	use debug_tools
	use pde_objs

	type(subdomain_str), intent(in) :: sub
	integer(kind=ikind), intent(in) :: itcount

	integer(kind=ikind) :: i, j, ii, jj, el, nd, xpos, xpos_loc
	real(kind=rkind) :: loc_error
	

	do i=1, ubound(sub%coarse_el,1)
	  j = sub%coarse_el(i)
	  loc_error = 0
	  do ii=1, ubound(ddcoarse_mesh(j)%elements,1)
	    el = ddcoarse_mesh(j)%elements(ii)
	    do jj = 1, ubound(elements%data,2)
	      nd = elements%data(el,jj)
	      xpos = pde(1)%permut(nd)
	      if (xpos > 0) then
		xpos_loc = sub%invpermut(xpos)
		loc_error = max(loc_error, abs(sub%xvect(xpos_loc,3)-sub%xvect(xpos_loc,2)) )
	      end if
	    end do
	  end do
	  if (loc_error > iter_criterion) then
	    ddcoarse_mesh(j)%iter_count = ddcoarse_mesh(j)%iter_count + 1
	  end if
	end do
  

      end subroutine search_error_cluster


      subroutine results_extractor()
	  use typy
	  use globals
	  use global_objs
	  use debug_tools
	  use pde_objs

	  integer(kind=ikind) :: i, j, proc


	  do proc=1, ubound(pde,1)
	      do i=1, nodes%kolik
		  if (pde(proc)%permut(i) > 0) then
		      pde(proc)%solution(i) = pde_common%xvect(pde(proc)%permut(i),3)
		  else
		      j = nodes%edge(i)
		      pde(proc)%solution(i) = pde(proc)%bc(j)%value
		  end if
	      end do
	  end do


	  pde_common%xvect(:,1) = pde_common%xvect(:,3)
	  pde_common%xvect(:,2) = pde_common%xvect(:,3)

      end subroutine results_extractor
 

      !> counts number of solved subdomains
      function domains_solved() result(no)
        use typy
        use globals
        use global_objs
        use decomp_vars
        
        integer(kind=ikind) :: no

        integer(kind=ikind) :: i

        no = 0

        do i=1, ubound(subdomain,1)
          if (abs(subdomain(i)%time - time+time_step) < 100*epsilon(time_step)) then
            no = no + 1
            subdomain(i)%solved = .true.
          end if
        end do

      end function domains_solved
      
            !> counts number of solved subdomains
      function ndofs_solved() result(no)
        use typy
        use globals
        use global_objs
        use decomp_vars
        
        integer(kind=ikind) :: no

        integer(kind=ikind) :: i

        no = 0

        do i=1, ubound(subdomain,1)
          if (subdomain(i)%solved) then
            no = no + subdomain(i)%ndof
          end if
        end do

      end function ndofs_solved
      
      subroutine correctdoms()
      
      end subroutine correctdoms
      
     subroutine getlocres(subdom, resvct)
      use typy
      use global_objs
      use pde_objs
      use globals
      use decomp_vars
      use debug_tools


      class(subdomain_str), intent(in) :: subdom
      real(kind=rkind), dimension(:), intent(out) :: resvct
      
      resvct = subdom%bvect - subdom%matrix%mul(subdom%xvect(:,3))
      
      
    
    
    end subroutine getlocres
      

    


end module subcycling
