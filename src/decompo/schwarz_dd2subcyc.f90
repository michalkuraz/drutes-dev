module schwarz_dd2subcyc
    use typy
    public :: schwarz_subcyc
    private :: build_xvect_loc
    private :: domains_solved
    private :: locmat_assembler
    private :: search_error_cluster
    private :: set_solved
    
    real(kind=rkind), dimension(:), allocatable, private, save :: resvct, corrvct
 
    contains

      !>solves the problem using Picard method with Schwarz multiplicative domain decomposition reflecting the nonlinear problem
      subroutine schwarz_subcyc(ierr, itcount, success)
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

	integer, intent(out) :: ierr
	integer(kind=ikind), intent(out) :: itcount
	integer(kind=ikind) ::  fin, pcg_it, subfin
	logical, intent(out) :: success
	integer(kind=ikind) :: i, proc,j, m, l, k, ii, nd
	real(kind=rkind) :: error, loc_error, res_error, reps, cumerr
	real(kind=rkind), dimension(:), allocatable, save :: coarse2fine
	logical :: dt_fine, first_run, picard_done, reset_domains
	integer :: ierr_loc
	integer(kind=ikind), dimension(:), allocatable :: radek, radek2, radek3
	character(len=128) :: text
	real(4), dimension(2,2) :: taaray
	real(4), dimension(2) :: rslt
	real(kind=rkind) :: val


        call etime(taaray(1,:), rslt(1))
	proc = ubound(pde,1)
	fin = maxval(pde(proc)%permut(:))
	itcount = 0
	inner_criterion = 1e-4

	!reset local-local cluster iteration count
        ddcoarse_mesh(:)%iter_count = 4
	
	!the residual vector is not allocated, thus we are at the beginning
	if (.not. allocated(resvct)) then
          first_run = .true.
	  allocate(resvct(fin))
	  allocate(corrvct(fin))          
          call set_subdomains()
        else
          first_run = .false.
        end if

        subdomain(:)%time_step = time_step
        
        subdomain(1)%time_step = time_step
        subdomain(1)%short_dt = .false.
        
        
        subdomain(:)%timeprev = time
        subdomain(:)%time = time 
        
              
                
        !reset local subdomain iteration count
        subdomain(:)%itcount = 0                

        !reset successes on domains
        subdomain(:)%solved = .false.
   	
	picard_done = .false.
	

        write(unit=terminal, fmt="(a)") "  "
        write(unit=terminal, fmt="(a, I4,a)") " Solving", ubound(subdomain,1),  " subdomains .... "
        
        schwarz: do

	  subdoms:  do i=1, ubound(subdomain,1)
	 
	    if (.not. subdomain(i)%solved) then    
	      call solve_subdomain(subdomain(i), reps=1e-10)
	    end if
     	    
	  end do subdoms
	 

	  
	  call wait()
	  
	  if (domains_solved() == ubound(subdomain,1)) then
	    call build_xvect_loc()
	    
	    do i=1, ubound(subdomain,1)
	      call combinevals(subdomain(i), short=.false.)
	    end do
	    if (domains_solved() == ubound(subdomain,1)) then
	      EXIT schwarz
	    end if
	  end if
	
	end do schwarz
	  
!   	    
!   	  call build_xvect()
!   
!   	    
! 	  call progressbar( int(100*ndofs_solved()/(1.0*ddinfo%ndofs_tot)))
! 	  
! 
! 	  if (domains_solved() == ubound(subdomain,1)) then
! 	    call etime(taaray(2,:), rslt(2))
! 	    do i=1, ubound(subdomain,1)
! 	      if (.not. subdomain(i)%critical .and. subdomain(i)%itcount>1 ) then
! 		reset_domains = .true.
! 		EXIT
! 	      else
! 		reset_domains = .false.
! 	      end if
! 	    end do
! 	   
! 	    call printtime("Time spent on iterations: ", rslt(2)-rslt(1))
! 	    
! 	    if (reset_domains) then
! 	      call build_xvect()
! 	      call set_subdomains()
! 	    end if
! 	    
! 	    
! 	    success = .TRUE.
! 	    ierr = 0
! 	    write(unit=file_ddinfo, fmt=*) time, itcount, "|", subdomain(:)%ndof
! 	    call flush(file_ddinfo)
! 	    EXIT picard_loop
! 	  end if
! 
! 	  if (itcount > max_itcount) then
! 	    call set_subdomains()
! 	    success = .FALSE.
! 	    ierr = 1
! 	    RETURN
! 	  end if
! 
!         end do picard_loop
! 
! 
! 
!         do proc=1, ubound(pde,1)
! 	  do i=1, ubound(subdomain,1)
! 	    dt_fine = pde(proc)%dt_check()
! 	  end do
!         end do
!         
!         if (.not.(dt_fine)) then
!           ierr = 2
!           success = .false.
!           RETURN
!         else
!           ierr = 0
!           call results_extractor()
!           do i=1, ubound(subdomain,1)
! 	    subdomain(i)%xvect(:,1) = subdomain(i)%xvect(:,3)
! 	  end do
!         end if

      end subroutine schwarz_subcyc
      
      subroutine solve_subdomain(sub, reps)
	use typy
	use decomp_vars
	use decomp_tools
	use sparsematrix
	use simplelinalg
	use pde_objs
	use debug_tools
	
	class(subdomain_str), intent(in out) :: sub
	real(kind=rkind), intent(in) :: reps
	
	integer(kind=ikind) :: subfin, i, s, nd, pnd, extdom, j, extpnd
	real(kind=rkind) :: error       
	integer :: ierr
	
        sub%itcount = 0


	sub%itcount = sub%itcount + 1


	call locmat_assembler(mtx=sub%matrix, bvect=sub%bvect, &
	      permut=sub%permut, dt=sub%time_step, invpermut=sub%invpermut, domain_id=sub%order, extended=.false., &
	      ierr=ierr)

	call locmat_assembler(mtx=sub%extmatrix, bvect=sub%extbvect, &
	      permut=sub%extpermut, dt=sub%time_step, invpermut=sub%extinvpermut, &
	      domain_id=sub%order, extended=.true.,ierr=ierr)


        call getres_loc(sub)

	resvct = 0
	
	subfin = sub%ndof
	      
	resvct(sub%permut(1:sub%ndof)) = sub%resvct%main
	      
	resvct(sub%extpermut(1:sub%extndof)) = sub%resvct%ext
	
	corrvct(1:subfin) = 0.0
	
	call diag_precond(a=sub%matrix, prmt=sub%permut,  mode=1)
			    
	call solve_matrix(sub%matrix, resvct, corrvct(1:subfin), ilev1=0, itmax1=subfin, reps1=reps)
	      
	call diag_precond(a=sub%matrix, x=corrvct(1:subfin), mode=-1)
		      
	error = maxval(abs(corrvct(1:subfin)))
	
	sub%xvect(:,2) = sub%xvect(:,2) + corrvct(1:subfin)
	  
	sub%xvect(:,3) = sub%xvect(:,2)
	
	do i=1, ubound(sub%xvect,1)
          pnd = sub%permut(i)
          nd = pde_common%invpermut(pnd)
          do j=1, ddinfo%nodesinextsub(nd)%pos
            extdom = ddinfo%nodesinextsub(nd)%data(j)
            extpnd = subdomain(extdom)%extinvpermut(pnd)
            if (extpnd /= 0) then
              subdomain(extdom)%extxvect(extpnd,2:3) = (subdomain(extdom)%extxvect(extpnd,2:3) + sub%xvect(i,2:3))/2.0_rkind
            end if
          end do
        end do
	
	
	if (error < iter_criterion) then
	  sub%time = sub%time + sub%time_step
	  if (abs(sub%time - time - time_step) < 100*epsilon(time)) then
	    call set_solved(sub, mode=1)
	  else
	    call set_solved(sub, mode=0)
	  end if
	end if
	

	
      
      end subroutine solve_subdomain

      
      subroutine combinevals(subdom, short)
	use typy
	use decomp_vars
	use pde_objs
	use debug_tools
	
	
	type(subdomain_str), intent(in out) :: subdom
	logical, intent(in) :: short
	integer(kind=ikind) :: i, glob, globgeom, j, sub, loc
	real(kind=rkind) :: value, oldvalue
	
	do i=1, ubound(subdom%xvect,1)
	  glob = subdom%permut(i)
	  globgeom = pde_common%invpermut(glob)
	  oldvalue = subdom%xvect(i,2)
	  if (ddinfo%nodesinsub(globgeom)%pos >= 1) then
	  
	    value = 0
	    do j=1, ddinfo%nodesinsub(globgeom)%pos
	      sub = ddinfo%nodesinsub(globgeom)%data(j)
	      loc = subdomain(sub)%invpermut(glob)
	      value = value + subdomain(sub)%returnval(loc, subdom%time)*prolong_mtx%get(glob,sub)
	    end do

	    if (short) then
	      if (subdom%short_dt) then
		subdom%xvect(i, 2:3) = value
	      end if
	    else
	      subdom%xvect(i, 2:3) = value
	    end if
	    
	    
	    if (subdom%solved .and. abs(value-oldvalue) > iter_criterion) then
	      call set_solved(subdom, mode=-1)
	    end if
	    
	  end if
	end do
	
      
      end subroutine combinevals
      
      
      subroutine build_xvect_loc()
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
	
	pde_common%xvect(:,1) = pde_common%xvect(:,3)
	
	do i=1, ubound(subdomain,1)
	  subdomain(i)%xvect(:,:) = pde_common%xvect(subdomain(i)%permut(1:subdomain(i)%ndof),:)
          subdomain(i)%extxvect(:,:) = pde_common%xvect(subdomain(i)%extpermut(1:subdomain(i)%extndof),:)
	end do
	
	
      end subroutine build_xvect_loc

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
          if (subdomain(i)%solved) then
            no = no + 1
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
      
    
   subroutine locmat_assembler(mtx, bvect, permut,  dt, invpermut, domain_id, extended, ierr)
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
      use decomp_vars
      
      class(extsmtx), intent(in out) :: mtx
      real(kind=rkind), dimension(:), intent(out) :: bvect
      integer(kind=ikind), dimension(:), intent(in) :: permut
      real(kind=rkind), intent(in) :: dt
      integer(kind=ikind), dimension(:), intent(in) :: invpermut
      integer(kind=ikind), intent(in) :: domain_id
      logical, intent(in) :: extended
      integer, intent(out) :: ierr
      integer(kind=ikind) :: el,j,k,l, proc, ll, limits, nd, ii, pnd, nnd, pnnd, coarseel, iii, i
      logical, dimension(:), allocatable, save :: elsolved
      type(integpnt_str) :: quadpnt
      real(kind=rkind) :: value

      
      if (.not. allocated(elsolved)) then
        allocate(elsolved(elements%kolik))
      end if
      
      
      call null_problem(mtx, bvect)

      
      elsolved = .false. 
      
!       limits = ubound(stiff_mat,1)/ubound(pde,1)
      
      proc = 1
      

      loop_nodes: do pnd=1, ubound(permut,1)

                      if (permut(pnd) == 0) then
                        EXIT loop_nodes
                      end if

                      nd = permut(pnd)

                      nd = pde_common%invpermut(nd)


                      do ii=1,nodes%el2integ(nd)%pos

                        el = nodes%el2integ(nd)%data(ii)

                        if (.not. elsolved(el)) then

                          quadpnt%type_pnt = "ndpt"
                          quadpnt%column = 1
                          quadpnt%ddlocal=.true.
                          coarseel = ddinfo%elincoarse(el)
        
                          
                          if (domain_id /= ddinfo%coarseinsub(coarseel)) then
                            quadpnt%extended = .false.
                            quadpnt%subdom = ddinfo%coarseinsub(coarseel)
                          else
                            quadpnt%subdom=domain_id
                            quadpnt%extended = extended
                          end if

                          do k = 1, ubound(elements%data,2)
                          
                            nnd = elements%data(el,k)
                            pnnd = pde(proc)%permut(nnd)
                            
                           
                            if (pnnd /= 0) then 
                              if (domain_id /= quadpnt%subdom) then                                   
                                quadpnt%globtime = .false.
                                quadpnt%time4eval = subdomain(domain_id)%time
                                quadpnt%column = 2
                              else
                                quadpnt%globtime = .true.
                              end if
                            end if
                            
                            quadpnt%order = nnd
                            elnode_prev(k) = pde(proc)%getval(quadpnt)
                          end do 
                          
			  
         
                          
                          call build_bvect(el, dt, quadpnt_in=quadpnt)

                          call build_stiff_np(el, dt, quadpnt_in=quadpnt)

                          call pde_common%time_integ(el, quadpnt_in=quadpnt)
                               
                          quadpnt%element = el
			  quadpnt%column = 2
			  quadpnt%type_pnt = "gqnd"
			  
			  quadpnt%ddlocal = .true.			  
			
			  stiff_mat = stiff_mat + cap_mat
			  
	        	  call in2global(el,mtx, bvect, invpermut)
 			  
                          elsolved(el) = .true.
                        end if
                      end do

      end do loop_nodes
     

    end subroutine 
    
    
    subroutine set_solved(sub, mode)
      use decomp_vars
      
      type(subdomain_str), intent(in out) :: sub
      !> 0 - the short time step is solved
      !! 1 - the long time step is solved
      !! 2 - it's done burning the bridges, no way back
      !! -1 - status solved is switched to unsolved
      integer, intent(in) :: mode
      
      select case(mode)
	case(0)
	  sub%xvect(1,:) = sub%xvect(3,:)
	case(1)
	  sub%xvect(1,:) = sub%xvect(3,:)
	  sub%solved = .true.
	case(2)
	  sub%xvect(1,:) = sub%xvect(3,:)
	  sub%xvect(2,:) = sub%xvect(3,:)
	  sub%xvect(4,:) = sub%xvect(3,:)
	  sub%solved = .true.
	case(-1)
	  sub%solved = .false.
	  sub%xvect(1,:) = sub%xvect(4,:)
	case default
	  print *, "code bug, incorrect mode value, exited from schwarz_dd2subcyc::set_solved"
	  error stop
      end select
    
    end subroutine set_solved


end module schwarz_dd2subcyc
