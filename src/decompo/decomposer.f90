module decomposer
  use typy
  public :: set_subdomains
  public :: init_decomp
  private :: make_quadmesh !, set_subdomains_old
  private :: read_coarse_mesh


  
  contains

    subroutine init_decomp()
      use typy
      use globals
      use global_objs
      use pde_objs
      use decomp_vars
      use geom_tools
      use core_tools
      use decomp_tools
      use printtools
      use readtools
      use debug_tools
    

      real(kind=rkind) :: coarse_density
      integer(kind=ikind) :: i,j,k, counter, l, m, n
      real(kind=rkind), dimension(2) :: bod, a, b, c
      real(kind=rkind), dimension(:,:), allocatable :: domain
      integer(kind=ikind), dimension(:), allocatable :: elcounter
      logical, dimension(:), allocatable :: cluster_def
      logical, dimension(:), allocatable :: node_counted
      logical :: exited
      integer :: ierr
      integer :: file_conf
      integer(kind=ikind) :: method

      call find_unit(file_conf,100)

      open(unit=file_conf, file="drutes.conf/mesh/coarse.conf", action="read", status="old", iostat=ierr)

      if (ierr /= 0) then
        print *, "unable to open file drutes.conf/mesh/coarse.conf, called from decomposer::init_decomp"
        error stop
      end if
      
      call fileread(ddinfo%overlaps, file_conf, ranges=(/0_ikind, huge(1_ikind)/))

      call fileread(i=method, fileid=file_conf, ranges=(/0_ikind,1_ikind/))
      
      call fileread(coarse_density, file_conf)


     select case(method)
        case(0)
            call make_quadmesh(coarse_density)
        case(1)
            call read_coarse_mesh()
      end select

      
      call fileread(use_coarselev, file_conf)
            
      allocate(domain(ubound(coarse_elements%data,2), drutes_config%dimen))

      allocate(elcounter(elements%kolik))

      allocate(elements%gc(elements%kolik, 2))

      allocate(cluster_def(elements%kolik))
      
      allocate(node_counted(nodes%kolik))

      i = ubound(pde,1)

      counter = maxval(pde(i)%permut(:))

      allocate(ddinfo%nd_dsjnt%data(counter))

      allocate(ddinfo%nd_dsjnt%crit(counter))
      
      allocate(elements%subdom(elements%kolik))

      cluster_def = .false.

      do i=1, elements%kolik
	a = nodes%data(elements%data(i,1),:)
	b = nodes%data(elements%data(i,2),:)
	c = nodes%data(elements%data(i,3),:)
	elements%gc(i, 1:drutes_config%dimen) = gravity_center(a,b,c)
      end do

      call write_log("filling the clusters and creating graph of the coarse mesh...")
	    
      do i=1, coarse_elements%kolik
	counter = 0
	do k=1, ubound(coarse_elements%data,2)
	  domain(k,:) = coarse_nodes%data(coarse_elements%data(i,k),:)
	end do
	call progressbar(int(100*i/coarse_elements%kolik))
	do j=1, elements%kolik
	  if (.not. cluster_def(j)) then
	    if (inside(domain, elements%gc(j,1:2))) then
	      counter = counter + 1
	      elcounter(counter) = j
	      cluster_def(j) = .true.
	    end if
	  end if
	end do

	if (counter > 0) then
	  if (counter < 4) call write_log("WARNING, your coarse mesh defines cluster with less then 4 fine mesh elements")
	  allocate(ddcoarse_mesh(i)%elements(counter))
	  ddcoarse_mesh(i)%elements = elcounter(1:counter)
	  ddcoarse_mesh(i)%el_count = counter
	else
	  call file_error(file_conf, message="The coarse mesh contains elements that does not cover the &
	  fine mesh elements, the COARSE mesh has INCORRECT geometry.")
	end if
      end do

      ddcoarse_mesh(:)%iter_count = 0
      ddcoarse_mesh(:)%subdomain = 0
      ddcoarse_mesh(:)%prevsub = 0
      
      do i=1, coarse_elements%kolik
        k = elements%material(ddcoarse_mesh(i)%elements(1),1)
        exited = .false.
	fine:  do j=1, ubound(ddcoarse_mesh(i)%elements,1)
	    l = ddcoarse_mesh(i)%elements(j)
	    if (k /= elements%material(l,1)) then
	      coarse_elements%material(i,1) = 0
	      exited = .true.
	      EXIT fine
	    end if
	end do fine
	 
	if (.not.exited) then
	  coarse_elements%material(i,1) = k
	end if
      end do
 

      call write_log("clusters filled.")
      
      deallocate(domain)
      deallocate(elcounter)
      deallocate(cluster_def)
      deallocate(node_counted)
            
      call find_unit(file_decomp)
      
      open(unit=file_decomp, file="out/dd/ddcount", action="write", status="replace", iostat=ierr)
      if (ierr /= 0) then
        ierr=system("mkdir out/dd")
        if (ierr /=0) then
          error stop "ERROR: unable to create directory out/dd, called from decomposer::init_decomp"
        else
          open(unit=file_decomp, file="out/dd/ddcount", action="write", status="replace", iostat=ierr)
          if (ierr /=0) then
            error stop "ERROR: unable to open file out/dd/ddcount, called from decomposer::init_decomp"
          end if
        end if
      end if
      
      call find_unit(file_ddinfo)
      
      open(unit=file_ddinfo, file="out/dd/ddinfo", action="write", status="replace")
      
      call init_addlevel()
      

        

    end subroutine init_decomp

 

    !> the subdomains are combined based on the gradients of the solution and the Picard iterations
    subroutine set_subdomains()
      use typy
      use globals
      use global_objs
      use pde_objs
      use decomp_vars
      use decomp_tools
      use geom_tools
      use core_tools
      use debug_tools
      use postpro


      logical :: remeshed
      integer(kind=ikind) :: i,j, k, l, counter, tmp_int, coarse_count, ii, subfin, processes_count, m, n, level_nei, decrease
      integer(kind=ikind) :: el, nd, sub
      integer(kind=ikind) :: level, currentlev, row
      real(kind=rkind)    :: avgval, tmp, tmp1, tmp2, dgdz, dgdx
      real(kind=rkind), dimension(2) :: locders
      real(kind=rkind), dimension(3,3) :: points
      logical :: exited
      integer(kind=ikind), dimension(:), allocatable, save :: subdomsdim, coarse_el
      integer(kind=ikind), dimension(:), allocatable, save :: subdoms
      logical, dimension(:), allocatable, save :: node_set
      integer(kind=ikind), dimension(:), allocatable, save :: domain_permut
      

      ddcoarse_mesh(:)%active = .true.
      ddcoarse_mesh(:)%evaluated = .false.
      

      processes_count = ubound(pde,1)
      if (.not.(allocated(subdomsdim))) then
	allocate(subdomsdim(ubound(ddcoarse_mesh,1)))
	allocate(subdoms(elements%kolik))
        allocate(coarse_el(coarse_elements%kolik))
        allocate(domain_permut(coarse_elements%kolik))
        allocate(node_set(maxval(pde(processes_count)%permut(:))))
      end if
      
      remeshed = .false.

      n=0


      !check if the cluster is critical (=steep derivatives or more than a single Picard iteration) 
      do i=1, coarse_elements%kolik
        row = 0
	exited = .false.
        incoarse: do j=1, ubound(ddcoarse_mesh(i)%elements,1)
                    l = ddcoarse_mesh(i)%elements(j)
		    row=row+1
		    do k=1,3
		      points(k,1:2) = nodes%data(elements%data(l,k),:)
		      points(k,3) = pde(1)%solution(elements%data(l,k))
		    end do
		    
		    call plane_derivative(points(1,:), points(2,:), points(3,:), dgdx, dgdz)
		    
!                     if ( ddcoarse_mesh(i)%iter_count > 1 .and. sqrt(dgdx*dgdx + dgdz*dgdz) > 0.0001 ) then
                      if ( ddcoarse_mesh(i)%iter_count > 1 ) then
		      ddcoarse_mesh(i)%treat_alone = .true.
		      exited = .true.
		      EXIT incoarse
		    end if
        end do incoarse
	  
	if (.not.exited) then
	  ddcoarse_mesh(i)%treat_alone = .false.
	end if
	if (coarse_elements%material(i,1) == 0) then
	  ddcoarse_mesh(i)%treat_alone = .true.
	end if
      end do

      
      currentlev = 1
      subdomsdim = 0
      ddcoarse_mesh(:)%prevsub = ddcoarse_mesh(:)%subdomain
      ddcoarse_mesh(:)%subdomain = 0


      do i=1, coarse_elements%kolik
        if (ddcoarse_mesh(i)%subdomain == 0) then
          ddcoarse_mesh(i)%subdomain = currentlev
          currentlev = currentlev + 1
        else if (ddcoarse_mesh(i)%subdomain < currentlev) then
          continue
        end if

        if (.not. ddcoarse_mesh(i)%treat_alone) then
	    do j=1, ubound(coarse_elements%neighbours,2)
              k = coarse_elements%neighbours(i,j) 
              if (k > 0) then 
		if (coarse_elements%material(i,1) == coarse_elements%material(k,1)) then
		  if (.not. ddcoarse_mesh(k)%treat_alone .and. ddcoarse_mesh(k)%subdomain == 0) then
		    ddcoarse_mesh(k)%subdomain= ddcoarse_mesh(i)%subdomain
		  end if
		  if (ddcoarse_mesh(k)%subdomain < ddcoarse_mesh(i)%subdomain .and. .not. ddcoarse_mesh(k)%treat_alone) then
		    level_nei = ddcoarse_mesh(k)%subdomain
		    decrease = ddcoarse_mesh(i)%subdomain
		    where (ddcoarse_mesh(:)%subdomain == decrease)
		      ddcoarse_mesh(:)%subdomain = level_nei
		    end where
		  end if
		end if
              end if
            end do
          end if
        end do
          
     
       currentlev = 0
       domain_permut = 0
       do i=1, ubound(ddcoarse_mesh,1)
        if (ddcoarse_mesh(i)%subdomain /= currentlev) then
          if (domain_permut( ddcoarse_mesh(i)%subdomain ) == 0) then
            currentlev = currentlev + 1
            domain_permut( ddcoarse_mesh(i)%subdomain ) = currentlev 
            ddcoarse_mesh(i)%subdomain = currentlev
          else
            ddcoarse_mesh(i)%subdomain = domain_permut( ddcoarse_mesh(i)%subdomain ) 
          end if
        end if
        if (ddcoarse_mesh(i)%subdomain /= ddcoarse_mesh(i)%prevsub ) then
          remeshed = .true.
        end if
      end do



      ddinfo%number = currentlev 
      ddinfo%t_change = time
      call cpu_time(ddinfo%cput_change)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!! begin if remeshed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ubound(subdomain,1) /= ddinfo%number .or. remeshed) then
	if (allocated(subdomain)) then
	  do i=1, ubound(subdomain,1) 
            call clear_subdomain(subdomain(i))
	  end do
	  deallocate(ddinfo%nodesinsub)
	  deallocate(ddinfo%nodesinextsub)
	  deallocate(subdomain)
	end if	
	
	do i=1, ubound(ddinfo%nd_dsjnt%crit,1)
	  call ddinfo%nd_dsjnt%crit(i)%clear()
	end do
	 
	allocate(subdomain(ddinfo%number))
        allocate(ddinfo%nodesinsub(nodes%kolik))
        allocate(ddinfo%nodesinextsub(nodes%kolik))
       


	do i=1, ubound(subdomain,1)
	  counter = 0
          coarse_count = 0
	  do j=1, coarse_elements%kolik
	    if (ddcoarse_mesh(j)%subdomain == i) then
	      tmp_int = counter + 1
	      counter = counter + ubound(ddcoarse_mesh(j)%elements,1)
              coarse_count = coarse_count + 1
              coarse_el(coarse_count) = j
	      subdoms(tmp_int:counter) = ddcoarse_mesh(j)%elements
	    end if
	  end do
          ! write fine elements
! 	  allocate(subdomain(i)%el(counter))
!  	  subdomain(i)%el = subdoms(1:counter)
	  do j=1, counter
	    call subdomain(i)%elsub%fill(subdoms(j))
	  end do
          !write coarse elements
          allocate(subdomain(i)%coarse_el(coarse_count))
          subdomain(i)%coarse_el(1:coarse_count) = coarse_el(1:coarse_count)
	end do

        ! create disjoint map of subdomain elements
	do i=1, ubound(subdomain,1)
	  do j=1, subdomain(i)%elsub%pos
	    k = subdomain(i)%elsub%data(j)
	    call elements%subdom(k)%fill(i)
	  end do
	end do

! 	call printmtx(elements%subdom) ; call wait()
        subdomain(:)%critical = .false.

        do i=1, coarse_elements%kolik
          if (ddcoarse_mesh(i)%treat_alone) then
            subdomain(ddcoarse_mesh(i)%subdomain)%critical = .true.
          end if
        end do
        
        
	! and now you can add desired number of levels
	call addlevel()

	
	! finish the subdomain initialization
	do i=1, ubound(subdomain,1)
	  call init_subdomain(subdomain(i), i)
	end do
	call write_log(text="domain decomposition was adapted to the new solution, number of subdomains is:",&
			int1=ddinfo%number, text2="@ simulation time:", real2=ddinfo%t_change)



       ! create a list of nodes with a list of subdomains, where the node belongs.
        do sub=1, ubound(subdomain,1)
	  do i=1, subdomain(sub)%elsub%pos
	    el = subdomain(sub)%elsub%data(i)
	    do j=1, ubound(elements%data,2)
	      nd = elements%data(el,j)
	      call ddinfo%nodesinsub(nd)%nrfill(sub)
	    end do
	  end do
	end do
	
	

      
	write(unit=file_decomp, fmt=*) time, ubound(subdomain,1)
	call flush(file_decomp)
	
	
	node_set = .false.
    
	!> create disjoint map of subdomains      
	do i=1, ubound(subdomain,1)
	  do j=1, subdomain(i)%elsub%pos
	    k = subdomain(i)%elsub%data(j)
	    do l=1, ubound(elements%data,2)
	      m = elements%data(k,l)
	      n = pde(1)%permut(m)
	      if (n > 0) then
		if (.not. node_set(n) .or. subdomain(i)%critical) then
		  ddinfo%nd_dsjnt%data(n) = i
		  node_set(n) = .true.
		  if (subdomain(i)%critical) then
		    call ddinfo%nd_dsjnt%crit(n)%fill(i)
		  end if
		end if
	      end if
	    end do
	  end do
	end do   
	
	do i=1, ubound(ddinfo%nd_dsjnt%data,1)
	  j=ddinfo%nd_dsjnt%data(i)
	  call subdomain(j)%disjoint%fill(i)
	end do
      end if

      ddinfo%ndofs_tot = 0
      do i=1, ubound(subdomain,1)
	ddinfo%ndofs_tot = ddinfo%ndofs_tot + subdomain(i)%ndof
      end do
      
      do i = 1, ubound(subdomain,1)
	subdomain(i)%order=i
      end do


    end subroutine set_subdomains
  
    
    subroutine init_subdomain(sub, id)
      use typy
      use sparsematrix
      use decomp_vars
      use globals
      use global_objs
      use pde_objs
      use core_tools
      use debug_tools

      
      type(subdomain_str), intent(in out) :: sub
      integer(kind=ikind), intent(in) :: id
      integer(kind=ikind) :: i,j,k,l, counter
      logical, dimension(:), allocatable, save :: wrklist
    
      if (.not.(allocated(wrklist))) then
	allocate(wrklist(maxval(pde(1)%permut(:))))
      end if
      
      allocate(sub%permut(3_ikind*sub%elsub%pos))
      allocate(sub%invpermut(maxval(pde(1)%permut(:))))
      
      wrklist = .false.
      sub%ndof = 0
      sub%permut = 0
      
      do i=1,sub%elsub%pos
	do j=1, ubound(elements%data,2)
	  l = sub%elsub%data(i)
	  k = pde(1)%permut(elements%data(l,j))
	  if (k /=0 ) then
	    wrklist(k) = .true.
	  end if
	end do
      end do
      
      
      counter = 0
      sub%invpermut = 0
      do i=1, ubound(wrklist,1)
	if (wrklist(i)) then
	  sub%ndof = sub%ndof + 1
	  sub%permut(sub%ndof) = i
          counter = counter + 1
          sub%invpermut(i) = counter
	end if
      end do
      

      call sub%matrix%init(maxval(pde(1)%permut(:)),sub%ndof)
      
      !!if subcycling      
      if (drutes_config%it_method == 2) then
	
	i = sub%ndextra%pos
		
	allocate(sub%extpermut(i))
	
	allocate(sub%extinvpermut(maxval(pde(1)%permut(:))))	
		
	sub%extinvpermut = 0
	
	sub%extpermut = 0
	
	j=1	
	do i=1, sub%ndextra%pos
	  if (pde(1)%permut(sub%ndextra%data(i)) > 0) then
	    sub%extpermut(j) =  pde(1)%permut(sub%ndextra%data(i))
	    j=j+1
	  end if	  
	end do
	
	sub%extndof = j-1
	
	
	call sub%extmatrix%init(maxval(pde(1)%permut(:)), sub%extndof)

	
	allocate(sub%extbvect(sub%extndof))
	
	l = 0
	do i=1, sub%ndextra%pos
	
 	  j = sub%ndextra%data(i)
 	  
 	  k = pde(1)%permut(j)

 	  if (k>0) then
 	    l = l + 1
	    sub%extinvpermut(k) = l
	  end if
	  
	end do
	
	
	allocate(sub%resvct%main(sub%ndof))
	
	allocate(sub%resvct%ext(sub%extndof))
	
	allocate(sub%extxvect(sub%extndof,4))
	
        sub%extxvect(:,1:3) = pde_common%xvect(sub%extpermut(1:sub%extndof),1:3)
        sub%extxvect(:,4) = pde_common%xvect(sub%extpermut(1:sub%extndof),1)

      end if
      
      allocate(sub%xvect(sub%ndof,4))
      

      
      allocate(sub%ovect(sub%ndof))
      
      allocate(sub%bvect(sub%ndof))

      
      sub%xvect(:,1:3) = pde_common%xvect(sub%permut(1:sub%ndof),1:3)
      
      sub%xvect(:,4) = pde_common%xvect(sub%permut(1:sub%ndof),1)
      

            
    end subroutine init_subdomain
    
    
    subroutine clear_subdomain(sub)
      use typy
      use sparsematrix
      use decomp_vars
      use core_tools
      use solver_interfaces
      
      type(subdomain_str), intent(in out) :: sub
      integer(kind=ikind) :: i
      
      call null_problem(sub%matrix)
      
      
      if (allocated(sub%xvect)) then
	deallocate(sub%xvect)
      else
        call write_log("strange, subdomain%xvect not allocated, called from decomposer::clear_subdomain")
      end if 
     
 
      call sub%elsub%clear(full = .true.)

      if (allocated(sub%coarse_el)) then
        deallocate(sub%coarse_el)
      else
        call write_log("strange, subdomain%coarse_el not allocated, called from decomposer::clear_subdomain")
      end if 
      
      if (allocated(sub%ovect)) then
	deallocate(sub%ovect)
      else
        call write_log("strange, subdomain%ovect not allocated, called from decomposer::clear_subdomain")
      end if 
      
      if (allocated(sub%permut)) then
	deallocate(sub%permut)
      else
        call write_log("strange, subdomain%permut not allocated, called from decomposer::clear_subdomain")
      end if 
      
      if (allocated(sub%invpermut)) then
	deallocate(sub%invpermut)
      else
        call write_log("strange, subdomain%invpermut not allocated, called from decomposer::clear_subdomain")
      end if
      
      if (allocated(sub%bvect)) then
        deallocate(sub%bvect)
      else
        call write_log("strange, subdomain%bvect not allocated, called from decomposer::clear_subdomain")
      end if

      call sub%disjoint%clear(full = .true.)

      do i=1, elements%kolik
	call elements%subdom(i)%clear()
      end do
      
      if (allocated(sub%extbvect)) then
	deallocate(sub%extbvect)
	deallocate(sub%extpermut)
	deallocate(sub%extinvpermut)
      end if
      
      if (allocated(sub%resvct%main)) then
	deallocate(sub%resvct%main)
      end if
      
      if (allocated(sub%resvct%ext)) then
	deallocate(sub%resvct%ext)
      end if
      
      
      call sub%elextra%clear(full = .true.)
      
      call sub%ndextra%clear(full = .true.)
      
      call null_problem(sub%extmatrix)
      
    
    end subroutine clear_subdomain


    subroutine read_coarse_mesh()
      use typy
      use global_objs
      use globals
      use decomp_vars
      use core_tools
      use readtools
      use debug_tools
      use decomp_tools
      use geom_tools

      integer :: fmesh, ierr, chi
      real :: ch
      integer(kind=ikind) :: i

      call find_unit(fmesh)
    
      open(unit=fmesh, file="drutes.conf/mesh/coarse.t3d", action="read", status="old", iostat=ierr)
  
      if (ierr /= 0) then
        print *, "unable to open file drutes.conf/mesh/coarse.t3d, called from decomposer::read_coarse"
        error stop
      end if

      call comment(fmesh)
      read(unit=fmesh, fmt=*) ch
      call comment(fmesh)
      read(unit=fmesh, fmt=*) coarse_nodes%kolik, ch, coarse_elements%kolik
      call comment(fmesh) 

      allocate(coarse_nodes%data(coarse_nodes%kolik, 2))
    
      allocate(coarse_elements%data(coarse_elements%kolik, 3))
      
      allocate(coarse_elements%material(coarse_elements%kolik,1))
      

      do i=1, coarse_nodes%kolik
        call comment(fmesh)
        read(unit=fmesh, fmt=*) ch, coarse_nodes%data(i,:)
      end do

      do i=1, coarse_elements%kolik
        call comment(fmesh)
        read(unit=fmesh, fmt=*) chi, coarse_elements%data(i,:)
      end do

      allocate(ddcoarse_mesh(coarse_elements%kolik))

      close(fmesh)

      call write_log("searching for neigbourhood clusters...")
      call find_neighbours(coarse_elements, coarse_nodes)

    end subroutine read_coarse_mesh
   
   
    !> subroutine that creates mesh of quadrangles -- coarse mesh (clusters)
    subroutine make_quadmesh(density)
      use typy
      use globals
      use global_objs
      use decomp_vars
      use printtools
      use core_tools
      use decomp_tools
      use geom_tools
      
      !> mesh density (approximatelly)
      real(kind=rkind), intent(in) :: density
      real(kind=rkind) :: xmax, xmin, zmax, zmin
      real(kind=rkind) :: xdens, zdens

      integer(kind=ikind) :: xcount, zcount, i,j, row

      xmax = maxval(nodes%data(:,1))
      xmin = minval(nodes%data(:,1))
      zmax = maxval(nodes%data(:,2))
      zmin = minval(nodes%data(:,2))
      if (((xmax - xmin)/density - int((xmax - xmin)/density)) > epsilon(xmax)) then
	xcount = int((xmax - xmin)/density) + 1
      else
	xcount = int((xmax - xmin)/density)
      end if

      xdens = (xmax - xmin)/(xcount*1.0_rkind)
 
      if (((zmax - zmin)/density - int((zmax - zmin)/density)) > epsilon(zmax)) then
	zcount = int((zmax - zmin)/density) + 1
      else
	zcount = int((zmax - zmin)/density)
      end if

      zdens = (zmax - zmin)/(zcount*1.0_rkind)
      
      allocate(coarse_nodes%data((xcount + 1)*(zcount + 1), 2))
      coarse_nodes%kolik = (xcount + 1)*(zcount + 1)

      allocate(coarse_elements%data(xcount*zcount, 4))
      coarse_elements%kolik = xcount*zcount

      allocate(ddcoarse_mesh(coarse_elements%kolik))
      
      allocate(coarse_elements%material(coarse_elements%kolik,1))
      
      allocate(coarse_elements%neighbours(coarse_elements%kolik, ubound(coarse_elements%data,2)))


      row = 0
      do j=0, zcount
	do i=0, xcount
	  row = row+1
	  coarse_nodes%data(row, 1) = xmin + i*xdens
	  coarse_nodes%data(row, 2) = zmin + j*zdens
	end do
      end do
   
      row = 0
      do j=1, zcount
	do i=1, xcount
	  row = row + 1
	  coarse_elements%data(row,1) = (xcount+1)*(j-1) + i
	  coarse_elements%data(row,2) = (xcount+1)*(j-1) + i + 1
	  coarse_elements%data(row,3) = (xcount+1)*j + i + 1
	  coarse_elements%data(row,4) = (xcount+1)*j + i
	end do
      end do

      
      print *, allocated(coarse_elements%neighbours) 
      
      call write_log("searching for neigbourhood clusters...")
      call find_neighbours(coarse_elements, coarse_nodes)
  

    end subroutine make_quadmesh

end module decomposer