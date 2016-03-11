! Copyright 2008 Michal Kuraz, Petr Mayer

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
!> solver for the advection dispersion equation

module postpro
  use typy
  
  public :: make_print
  public :: write_obs
  public :: get_RAM_use
  private :: print_scilab, print_pure, print_gmsh
  


  contains
 
    subroutine make_print(behaviour, dt, name)
      use typy
      use globals
      use global_objs
      use core_tools
      use geom_tools
      use pde_objs
      use debug_tools

      character(len=*), intent(in), optional                :: behaviour
      real(kind=rkind), intent(in), optional                :: dt
      character(len=*), intent(in), optional                :: name
      logical                                               :: anime
      integer(kind=ikind)                                   :: mode, no_prints
      character(len=256), dimension(:,:), allocatable       :: filenames
      integer(kind=ikind)                                   :: i, i_err, j, layer, proc, run
      character(len=64)                                     :: forma
      integer, dimension(:,:), pointer, save                :: ids
      integer, dimension(:,:), allocatable, target, save    :: ids_obs
      integer, dimension(:,:), allocatable, target, save    :: ids_anime
      character(len=5)                                      :: extension
      character(len=15)                                     :: prefix
      real(kind=rkind)                                      :: distance, avgval, val1, val2, val3, tmp, flux
      logical                                               :: first_run
      type(integpnt_str), dimension(:), allocatable, save   :: quadpnt
      integer(kind=ikind), save                             :: anime_run, anime_dec
      integer                                               :: ierr
      character(len=512) :: filename, iaction


      if (present(name)) then
	select case(name)
	  case("obs")
	    anime = .false.
	  case("avi")
	    anime = .true.
	  case default
	    print *, "this is a BUG, contact developer, called from postpro::make_print()"
            error STOP
        end select
      else
	anime = .false.
      end if

      if (present(behaviour) ) then
	select case(behaviour)
	  case("all_in_one")
	    mode = -1
	  case("separately")
	    mode = 0
	  case default
	    print *, "Incorrect behaviour specified, called from postpro::make_print"
	    ERROR stop
	end select
      else
	mode = 0
      end if

      if (drutes_config%dimen < 2  .or. www) then
	extension = ".dat"
      else if (drutes_config%mesh_type /=3 ) then
	extension = ".sci"
      else
	extension = ".gmsh"
      end if

      call write_log("making output files for observation time")
      
      
      if (.not. allocated(quadpnt)) then
	select case(drutes_config%dimen)
	  case(1)
	    allocate(quadpnt(2))
	  case(2)
	    allocate(quadpnt(1))
	end select
	first_run = .true.
      else
	first_run = .false.
      end if

      
      if (.not. anime) then
	postpro_run = postpro_run + 1
      else
        if (first_run) then
	  anime_run = 0
	end if
	anime_run = anime_run + 1
      end if

      
      if (.not. anime) then
	do 
	  if (postpro_run/(10**postpro_dec) < 1) then
	    EXIT
	  else
	    postpro_dec = postpro_dec + 1
	  end if
	end do
      end if
  
      if (.not. allocated(ids_obs)) then
	allocate(ids(ubound(pde,1), 4))
	allocate(ids_obs(ubound(pde,1), 4))
	allocate(ids_anime(ubound(pde,1), 4))
      end if
      
      if (anime) then
	ids => ids_anime
      else
	ids => ids_obs
      end if
	

      allocate(filenames(ubound(pde,1),4))
 
      
      if (anime) then
	prefix = "out/anime/"
	write(unit=forma, fmt="(a, I7, a)") "(a, a, a, a, a, I", 1," a)"
	run = anime_run
      else
	prefix = "out/"
	write(unit=forma, fmt="(a, I7, a)") "(a, a, a, a, a, I", postpro_dec," a)"
	run = postpro_run
      end if

      do proc=1, ubound(pde,1)

	write(unit=filenames(proc,1), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", &
	trim(pde(proc)%solution_name(1)), "-",  run, trim(extension)


	write(unit=filenames(proc,2), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", trim(pde(proc)%solution_name(1)), &
						  "-el_avg-",  run,  trim(extension)

	write(unit=filenames(proc,3), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_",  trim(pde(proc)%mass_name(1)), "-", & 
					     run,  trim(extension)

	write(unit=filenames(proc,4), fmt=forma) trim(prefix), trim(pde(proc)%problem_name(1)), "_", trim(pde(proc)%flux_name(1)), "-", &
						  run,  trim(extension)
! 	
	if ( .not. anime .and. (mode == 0 .or. postpro_run == 1 ) .or. &
	  (anime_run == 1 .and. anime)) then

	  do i=1, ubound(filenames,2)
	    call find_unit(ids(proc, i), 6000)
	    open(unit=ids(proc, i), file=trim(filenames(proc,i)), action="write", status="replace", iostat=ierr)
	    if (ierr /= 0) then
	      call system("mkdir out/anime")
	      open(unit=ids(proc, i), file=trim(filenames(proc,i)), action="write", status="replace", iostat=ierr)
	      if (ierr /= 0) then
		print *, "unexpected system error, called from postpro::make_print()"
		error stop
	      end if
	    end if
	  end do  
	end if


	quadpnt(:)%type_pnt = "ndpt"
	
	select case(observe_info%method)
	  case(1)
	      quadpnt(:)%column = 3
	  case(2)
	      if (.not. present(dt) .and. time > epsilon(time) .and. .not. (abs(time-end_time) < epsilon(time) .or. time > end_time)) then
! 	       print *, abs(time-end_time), time, end_time
		print *, "this is a BUG, insufficient input parameters, contact developer, called from postpro::make_print()"
		error STOP
	      end if
	      if (time > epsilon(time) .and. abs(time-end_time) > epsilon(time) .and. time<end_time) then
		quadpnt(:)%column = 4
		pde_common%xvect(:,4) = (pde_common%xvect(:,3)-pde_common%xvect(:,1))/time_step*dt + pde_common%xvect(:,1)
	      else if (time < epsilon(time) ) then
		quadpnt(:)%column = 1
	      else if (abs(time-end_time) < epsilon(time) .or. time > end_time) then
		quadpnt(:)%column = 3
	      end if
		
	 end select
	 


	if (drutes_config%dimen == 1 .or. www) then
  
	  call print_pure(ids(proc,:), proc, quadpnt)
	  
	else if  (drutes_config%mesh_type <= 2 ) then

	  if (present(dt)) then
	    call print_scilab(ids(proc,:), proc, quadpnt, dt)
	  else
	    call print_scilab(ids(proc,:), proc, quadpnt)
	  end if
	  

	else

	  call print_gmsh(ids(proc,:), proc, quadpnt)
	
	end if
      end do


      if (.not. anime) then
   	  do proc=1, ubound(pde,1)
	    do i=1,ubound(ids,2)
	      close(ids(proc,i))
	    end do
	  end do	
      else      
	do proc=1, ubound(pde,1)
	  do i=1,ubound(ids,2)
	    write(unit=ids(proc,i), fmt="(a,I6.6,a)" ) "xs2png(0, 'K-", anime_run , ".png');"
	    write(unit=ids(proc,i), fmt=*) "clear"
	    write(unit=ids(proc,i), fmt=*) "   "
	    call flush(ids(proc,i))
	  end do
	end do
	
      end if
      
     if (anime) then
       inquire(unit=ids(1,1), name=filename, action=iaction) 
     end if

      deallocate(filenames)
    end subroutine make_print




    subroutine write_obs()
      use typy
      use globals
      use global_objs
      use geom_tools
      use pde_objs
      use debug_tools
		
      integer(kind=ikind)               :: i,j,m,l, k, n, layer, proc, D
      real(kind=rkind)	              :: val, avgval, massval
      real(kind=rkind), dimension(3)    :: advectval, gradient, advectval2
      real(kind=rkind), dimension(3)    :: a,b,c
      real(kind=rkind) 		      :: tmp, zdist, xdist, tmp1, tmp2, pointdist, dgdx, dgdz, dgdy
      integer(kind=ikind), dimension(3) :: nodes_id
      type(integpnt_str) :: quadpnt
      
      quadpnt%type_pnt = "obpt"
      quadpnt%column=3
      D = drutes_config%dimen
      

      do proc=1, ubound(pde,1)
	do i =1, ubound(observation_array,1)
	  layer = elements%material(observation_array(i)%element, proc)

	  quadpnt%order = i
	  
	  call pde(proc)%flux(layer=layer, quadpnt=quadpnt,  vector_out=advectval(1:D))

	  val = pde(proc)%getval(quadpnt)
	  
	  observation_array(i)%cumflux(proc) = observation_array(i)%cumflux(proc) + &
		    sqrt(dot_product(advectval(1:D), advectval(1:D)))*time_step
	  
	  write(unit=pde(proc)%obspt_unit(i), fmt="(50E24.12E3)") time, val,  advectval(1:D), observation_array(i)%cumflux(proc)

	  call flush(pde(proc)%obspt_unit(i))
	  
	end do
      end do
      
    end subroutine write_obs
 

  !> this subroutine parses /proc/[PID]/status file in order to get RAM consumption statistics, it works only under POSIX systems
  subroutine get_RAM_use()
    use typy
    use globals
    use core_tools
    integer :: PID, fileid, i, ierr
    integer(kind=ikind) :: bytes
    character(len=2) :: byte_unit
    character(len=7) :: ch
    character(len=256) :: filename, format

    PID = getpid()

    call find_unit(fileid)


    i = 0
    do
      i = i + 1
      if ((1.0*PID)/(10**i) < 1) then
	EXIT
      end if
    end do

    write(unit=format, fmt="(a, I7, a)") "(a, I", i, ", a)"

!     write(unit=format, fmt = *) "(I.", i, ")"

    write(unit=filename, fmt=format) "/proc/", PID, "/status"
    
    open(unit=fileid, file=filename, action="read", status="old", iostat=ierr)

    if (ierr /= 0) then
      write(unit=terminal) "WARNING! this is not POSIX system, unable to get RAM consumption"
      RETURN
    end if

    do 
      read(unit=fileid, fmt=*, iostat=ierr) ch
      if (ch == "VmPeak:") then
	backspace fileid
	EXIT
      end if

      if (ierr /=0) then
	print *, "unable to fetch memory consumption from system files"
	RETURN
      end if
    end do

    read(unit=fileid, fmt=*) ch, bytes, byte_unit

    call write_log(text="Peak RAM  consumption on image", int1=1_ikind*THIS_IMAGE(), text2="was:", int2=bytes, text3=byte_unit)

    do 
      read(unit=fileid, fmt=*) ch
      if (ch == "VmSwap:") then
	backspace fileid
	EXIT
      end if

      if (ierr /=0) then
	print *, "unable to fetch swap consumption from system files"
	RETURN
      end if
    end do

    read(unit=fileid, fmt=*) ch, bytes, byte_unit

    call write_log(text="Peak SWAP consumption on image", int1=1_ikind*THIS_IMAGE(), text2="was:", int2=bytes, text3=byte_unit)


    close(fileid)

  end subroutine get_RAM_use
  
  subroutine print_scilab(ids, proc, quadpnt, dt)
    use typy
    use globals
    use global_objs
    use pde_objs
    use geom_tools
    
    integer, dimension(:), intent(in) :: ids
    integer(kind=ikind), intent(in) :: proc
    type(integpnt_str), dimension(:), intent(in out) :: quadpnt
    real(kind=rkind), intent(in), optional :: dt
    integer(kind=ikind) :: i, j, layer
    real(kind=rkind) :: tmp, flux
    real(kind=rkind), dimension(3,8) :: body
    real(kind=rkind), dimension(2) :: vct1, vct2
    real(kind=rkind), dimension(8) :: vct_tmp
    real(kind=rkind), dimension(2,2) :: blabla
    real(kind=rkind), dimension(:), allocatable, save :: deriv
    integer(kind=ikind) :: time_dec
    real(kind=rkind) :: curtime
    
    
    if (present(dt)) then
      curtime = time - time_step + dt
    else
      curtime = time
    end if
  
    do i=1, ubound(ids,1)
      write(unit=ids(i), fmt=*) "//", curtime
      write(unit=ids(i), fmt=*) "nt =", elements%kolik, ";"
      write(unit=ids(i), fmt=*) "x=zeros(nt,3);"
      write(unit=ids(i), fmt=*) "y=zeros(nt,3);"
      write(unit=ids(i), fmt=*) "z=zeros(nt,3);"
    end do

    if (.not. allocated(deriv)) then
      allocate(deriv(drutes_config%dimen))
    end if


    do i=1, elements%kolik
      do j=1,ubound(elements%data,2)
	body(j,1:2) = nodes%data(elements%data(i,j),:)
	quadpnt(1)%order = elements%data(i,j)
	body(j,3) = pde(proc)%getval(quadpnt(1))
      end do
    
      ! mass (constant over element)
      layer = elements%material(i, proc)
      tmp = sum(body(:,3))/ubound(elements%data,2)
      call pde(proc)%pde_fnc(proc)%dispersion(pde(proc), layer, x=(/tmp/), tensor=blabla)
      body(:,4) = maxval(blabla)
      body(:,5) = pde(proc)%mass(layer, x=(/tmp/))

      call plane_derivative(body(1,1:3), body(2,1:3), body(3,1:3), deriv(1), deriv(2))

      call pde(proc)%flux(layer, x=(/tmp/), vector_in=deriv, scalar=flux)

      body(:,6) = flux

      
      vct1 = body(3, 1:2) - body(1, 1:2) 
      vct2 = body(3, 1:2) - body(2, 1:2)
      
      if ( (vct1(1)*vct2(2)-vct1(2)*vct2(1)) > 0.0_rkind) then
	vct_tmp = body(3,:)
	body(3,:) = body(2,:)
	body(2,:) = vct_tmp
      else
	CONTINUE
      end if

      do j=1, ubound(ids,1)
	write(unit=ids(j), fmt=*) "x(", i, ",1) =", body(1,1), ";"
	write(unit=ids(j), fmt=*) "x(", i, ",2) =", body(2,1), ";"
	write(unit=ids(j), fmt=*) "x(", i, ",3) =", body(3,1), ";"
	
	write(unit=ids(j), fmt=*) "y(", i, ",1) =", body(1,2), ";"
	write(unit=ids(j), fmt=*) "y(", i, ",2) =", body(2,2), ";"
	write(unit=ids(j), fmt=*) "y(", i, ",3) =", body(3,2), ";"

	write(unit=ids(j), fmt=*) "z(", i, ",1) =", body(1,2+j), ";"
	write(unit=ids(j), fmt=*) "z(", i, ",2) =", body(2,2+j), ";"
	write(unit=ids(j), fmt=*) "z(", i, ",3) =", body(3,2+j), ";"
      end do

    end do
  
    time_dec = 0
    if (curtime>epsilon(curtime)) then
      do
	if (curtime*10.0_rkind**time_dec > 1) then
	  EXIT 
	else
	  time_dec = time_dec + 1
	end if
      end do
    end if
    


    do i=1, ubound(ids,1)
      write(unit=ids(i), fmt=*) "f=gcf();"
      write(unit=ids(i), fmt=*) "clf(f,'reset');"
      write(unit=ids(i), fmt=*) "f.color_map=jetcolormap(256);"
      write(unit=ids(i), fmt=*) "colorbar(min(z),max(z));"
      if (time_dec < 2) then
	write(unit=ids(i), fmt="(a,F10.2,a,a,a,a)")  "xtitle('$\mathbf{\LARGE t= ", curtime,"  ",&
	      "} \quad \mbox{\Large ",   trim(time_units), "}$')"
      else
							    
!                                                                  a	                 F10.2        a
       write(unit=ids(i), fmt="(a,F10.2,a,a,I16,a,a,a)")  "xtitle('$\mathbf{\LARGE t= ", curtime*10**time_dec,"  ",&
!                    a            I16         a                         a                a
	      "\times 10^{-", time_dec,"}} \quad \mbox{\Large ",   trim(time_units), "}$')"
      end if
      
      

!   xtitle('$\mathbf{\LARGE t=     300.00  } \quad \mbox{\Large  hrs}$')
!       write(unit=ids(proc,i), fmt="(a,F10.2,a,a)")  "xtitle('$\mathbf{\LARGE t= ", time,"  ", "} \quad \mbox{\Large  hrs}$')"
      write(unit=ids(i), fmt=*) "plot3d1(x',y',z',alpha=0, theta=-90);"

    end do 
  
  end subroutine print_scilab
  
  
  subroutine print_pure(ids, proc, quadpnt)
    use typy
    use globals
    use global_objs
    use pde_objs
    
    integer, dimension(:), intent(in) :: ids
    integer(kind=ikind), intent(in) :: proc
    type(integpnt_str), dimension(:), intent(in out) :: quadpnt
    integer(kind=ikind) :: i, layer
    real(kind=rkind) ::  distance, flux, avgval

  
    do i=1, nodes%kolik
      quadpnt(1)%order = i
      write(unit=ids(1), fmt="(50E24.12E3)") nodes%data(i,:), pde(proc)%getval(quadpnt(1)) 
      layer = elements%material(nodes%element(i)%data(1), proc)

      call pde(proc)%flux(layer, quadpnt(1), scalar=flux)

      write(unit=ids(4), fmt=*) nodes%data(i,:), flux
    end do


    do i=1, elements%kolik
      distance = (nodes%data(elements%data(i,1),1) + nodes%data(elements%data(i,2),1))/2.0
      quadpnt(1)%order = elements%data(i,1)
      quadpnt(2)%order = elements%data(i,2)
      avgval = (pde(proc)%getval(quadpnt(1)) + pde(proc)%getval(quadpnt(2)))/2.0
      layer = elements%material(i, proc)
      write(unit=ids(2), fmt="(50E24.12E3)") distance, avgval
      avgval = (pde(proc)%getval(quadpnt(1)) + pde(proc)%getval(quadpnt(2)))/2.0
      write(unit=ids(3), fmt="(50E24.12E3)") distance, avgval

    end do
  
  end subroutine print_pure
  
  subroutine print_gmsh(ids, proc, quadpnt)
    use typy
    use globals
    use global_objs
    use pde_objs
    
    integer, dimension(:), intent(in) :: ids
    integer(kind=ikind), intent(in) :: proc
    type(integpnt_str), dimension(:), intent(in out) :: quadpnt
    integer(kind=ikind) :: i
  

    write(unit=ids(1), fmt=*) "$MeshFormat"
    write(unit=ids(1), fmt=*) "2.2 0 8"
    write(unit=ids(1), fmt=*) "$EndMeshFormat"
    write(unit=ids(1), fmt=*) "$Nodes"
    write(unit=ids(1), fmt=*) nodes%kolik
    
    do i=1, nodes%kolik
      write(unit=ids(1), fmt=*) i,  nodes%data(i,:)
    end do
    write(unit=ids(1), fmt=*) "$EndNodes"
    write(unit=ids(1), fmt=*) "$Elements"
    write(unit=ids(1), fmt=*) elements%kolik
      do i=1, elements%kolik
      write(unit=ids(1), fmt=*) i,  elements%data(i,:)
    end do
    write(unit=ids(1), fmt=*) "$EndElements"
    write(unit=ids(1), fmt=*) "$NodeData"
    write(unit=ids(1), fmt=*) nodes%kolik
    do i=1, nodes%kolik
      write(unit=ids(1), fmt=*) i, pde(proc)%solution(i)
    end do
    write(unit=ids(1), fmt=*) "$EndNodeData"
  
  end subroutine print_gmsh

  
    
end module postpro