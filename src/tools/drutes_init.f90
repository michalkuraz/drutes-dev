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

module drutes_init
  public  :: parse_globals, init_observe, init_measured, get_cmd_options
  private :: set_global_vars,  open_obs_unit, pde_constructor, init_obstimes


  contains

    !> routine to open all files to read and write
    subroutine parse_globals()
      use typy
      use globals
      use pde_objs
      use core_tools
      use read_inputs
      use readtools


      integer :: i_err, i
      character(len=256) ::  writer
      character(len=8192) :: dirname
      real(kind=rkind) :: const


      call find_unit(logfile, 2000)

      open(unit=logfile, file="out/DRUtES.log", action="write", status="replace", iostat=i_err)
      
      if (i_err /= 0) then
	i_err=system("mkdir out")
        open(unit=logfile, file="out/DRUtES.log", action="write", status="replace", iostat=i_err)
	if (i_err /= 0) then
	  print *, "unable to open and create directory out (is it UN*X type os?), called  from drutes_init::parse_globals"
	  ERROR STOP
	end if
	  
      end if
      
      write(unit=logfile, fmt=*)  "  _____  _____  _    _ _   ______  _____ "
      write(unit=logfile, fmt=*)  " |  __ \|  __ \| |  | | | |  ____|/ ____|"
      write(unit=logfile, fmt=*)  " | |  | | |__) | |  | | |_| |__  | (___  "
      write(unit=logfile, fmt=*)  " | |  | |  _  /| |  | | __|  __|  \___ \ "
      write(unit=logfile, fmt=*)  " | |__| | | \ \| |__| | |_| |____ ____) | "
      write(unit=logfile, fmt=*)  " |_____/|_|  \_\\____/ \__|______|_____/ "
      
      write(unit=logfile, fmt=*)  " "
      write(unit=logfile, fmt=*)  " "
      write(unit=logfile, fmt=*)  " "
    
      write(unit=logfile, fmt=*, iostat = i_err) "DRUtES is initialising, reading input files"

      call find_unit(file_itcg, 2000)

      open(unit=file_itcg, file="out/cgit.count", action="write", status="replace", iostat=i_err)
      
      call write_log(text="total number of parallel images is:", int1=1_ikind*NUM_IMAGES())

      call write_log("reading input files and allocating...")

      !global.conf
      call find_unit(file_global, 200)
      
      open(unit=file_global,file="drutes.conf/global.conf", action="read", status="old", iostat = i_err)


      if (i_err /= 0) then
	print *, "missing drutes.conf/global.conf file"
	ERROR STOP
      end if
      
      
      if (www) then
	if (.not. dirglob%valid) then
	  print *, "if you have specified option www you MUST specify location with global settings"
	  print *, "specify option --dir-global path/to/your/global/settings"
	  ERROR STOP
	end if
	call find_unit(file_wwwglob, 200)
	write(dirname, fmt=*) trim(dirglob%dir), "drutes_global.conf/global.conf" 
	open(unit=file_wwwglob,file=trim(adjustl(dirname)), action="read", status="old", iostat = i_err)
	if (i_err /= 0) then
	  print *, "incorrect definition of global confispozditguration file"
	  print *, "your definition was: ", trim(dirglob%dir)
	  print *, "the constructed path was: ", trim(adjustl(dirname))
	  ERROR stop
	end if

      end if
     
      
      call read_global()
     
      call set_global_vars()
      
      call init_obstimes()

       !read mesh data
      call find_unit(file_mesh, 200)
      select case(drutes_config%dimen)
        case(1)
          open(unit=file_mesh, file="drutes.conf/mesh/drumesh1d.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/mesh/drumesh1d.conf file"
            ERROR STOP
          end if
          call read_1dmesh_int()

        case(2)
	  select case(drutes_config%mesh_type)
	    case(1)
	      open(unit=file_mesh, file="drutes.conf/mesh/drumesh2d.conf", action="read", status="old", iostat=i_err)
	      if (i_err /= 0) then
		print *, "missing drutes.conf/mesh/drumesh2d.conf file"
		ERROR STOP
	      end if
	      call read_2dmesh_int()
	    case(2)
	      open(unit=file_mesh, file="drutes.conf/mesh/mesh.t3d", action="read", status="old", iostat=i_err)
	      if (i_err /= 0) then
		print *, "missing drutes.conf/mesh/mesh.t3d file"
		ERROR STOP
	      end if
	      call read_2dmesh_t3d()
	    case(3)
	      open(unit=file_mesh, file="drutes.conf/mesh/mesh.gmsh", action="read", status="old", iostat=i_err)
	      if (i_err /= 0) then
		print *, "missing drutes.conf/mesh/mesh.gmsh file"
		ERROR STOP
	      end if
	      call read_2dmesh_gmsh()
	    case default
	      write(unit=terminal, fmt=*) "INCORRECT mesh type definition"
	      write(unit=terminal, fmt=*)"you have defined value:", drutes_config%mesh_type
	      write(unit=terminal, fmt=*) "the available options are: 1 - internal mesh generator" 
	      write(unit=terminal, fmt=*) "   (very simple uniform meshes, for debuging only"
	      write(unit=terminal, fmt=*)"                           2 - t3d mesh generator"
	      write(unit=terminal, fmt=*)"                           3 - gmsh mesh generator"
	      write(unit=terminal, fmt=*)"exiting, called from drutes_init::parse_globals"
	      error STOP
	   end select
        case default
          write(unit=terminal, fmt=*)"ERROR: unsupported problem dimension, the specified dimension was: ", drutes_config%dimen
          write(unit=terminal, fmt=*)"currently only 1D and 2D is supported"
          ERROR STOP 
      end select

      call pde_constructor(pde_common%processes)
      
      if (debugmode) then
	call find_unit(debugfile)
	open(unit=debugfile, file="out/debug.infos", action="write", status="replace")
      end if

      

    end subroutine parse_globals


    subroutine pde_constructor(processes)
      use typy
      use global_objs
      use globals
      use pde_objs
      use dummy_procs

      integer(kind=ikind), intent(in) :: processes
      
      integer :: i, j

      allocate(pde(processes))

      do i=1, processes
	allocate(pde(i)%pde_fnc(processes))
	do j=1, processes
	  pde(i)%pde_fnc(j)%dispersion => dummy_tensor
	  pde(i)%pde_fnc(j)%convection => dummy_vector
	  pde(i)%pde_fnc(j)%der_convect => dummy_vector
	  pde(i)%pde_fnc(j)%reaction => dummy_scalar
	  pde(i)%pde_fnc(j)%zerord => dummy_scalar
	  pde(i)%pde_fnc(j)%elasticity => dummy_scalar
	end do  
        allocate(pde(i)%solution(nodes%kolik))
        allocate(pde(i)%obspt_unit(ubound(observation_array,1)))
        allocate(pde(i)%permut(nodes%kolik))
        pde(i)%mass => dummy_scalar
        pde(i)%flux => dummy_vector
        pde(i)%dt_check => time_check_ok
        pde(i)%process_change => do_nothing
        pde(i)%getval => getvalp1
        pde(i)%order = i
      end do


    end subroutine pde_constructor
    
    subroutine get_cmd_options()
      use globals
      use core_tools
      
      character(len=64) :: arg, oarg
      character(len=4096) :: dir
      integer :: i, fileid, ierr, j
      integer,  dimension(:), allocatable :: skip 
      
      allocate(skip(1))
      skip = 0
      
      do i = 1, iargc()
	call getarg(i, arg)
      
        if (maxval(skip) == 0) then
	  select case(trim(arg))
	    case("-o")
		    call getarg(i+1,oarg)
		    if (trim(oarg) == "www") then
		      www = .true.
		    end if
		    skip = 1
		    
	    case("--print-level")
	            call getarg(i+1,oarg)
		    if (trim(oarg) == "1") then
		      call find_unit(terminal, 2000)
		      open(unit=terminal, file="out/screen.log", action="write", status="replace", iostat=ierr)
		      if (ierr /=0) then
			call system("mkdir out")
			open(unit=terminal, file="out/screen.log", action="write", status="replace", iostat=ierr)
		      end if
		      if (ierr /=0) then
			print *, "ERROR: directory out cannot be created"
			ERROR STOP
		      end if
		    else if (trim(oarg) == "0") then
		      CONTINUE
		    else if (trim(oarg) == "2") then
			call find_unit(terminal, 2000)
			open(unit=terminal, file="/dev/null", action="write", status="replace", iostat=ierr)
		        if (ierr /=0) then
			  print *, "ERROR: file /dev/null does not exist, are you running GNU/Linux based os?"
			  print *, "       if not do not set --print-level 2 in your command line option"
			  ERROR STOP
			end if
		      else
			print *, "incorrect argument after --print-level, exiting..."
			ERROR STOP
		      end if
		    skip = 1
		    terminal_assigned = .true.
				    
	    case("-v")
		print *, " " //achar(27)//'[94m', "DRUtES" //achar(27)//'[0m', &
		" version: " //achar(27)//'[92m', version_id, " " //achar(27)//'[0m'
    ! 		print *, "compiled with:", compiler_version()
    ! 		print *, "with options:", compiler_options()
		    STOP
	    case("--tmax")
		  call getarg(i+1, oarg)
		  read(unit=oarg, fmt=*, iostat=ierr) cpu_max_time
		  if (ierr /= 0 .or. cpu_max_time <= 0) then
		   print *, "incorrect maximal cpu time value (the value after --tmax option)"
		   ERROR STOP
		  end if
		  call getarg(i+2, oarg)
		  select case(trim(oarg))
		    case("s")
		      CONTINUE
		    case("min")
		      cpu_max_time = cpu_max_time * 60
		    case("hrs")
		      cpu_max_time = cpu_max_time * 3600
		    case("day")
		      cpu_max_time = cpu_max_time * 86400
		    case default
		      print *, "you have specified incorrect time units for maximal CPU time (--tmax option)"
		      print *, "the available units are [s, min, hrs, day]"
		      ERROR STOP
		  end select
		  cpu_time_limit = .true.
		  deallocate(skip)
		  allocate(skip(2))
		  skip = 1
	    case("--dir-local")
		  call getarg(i+1, dir)
		  call chdir(trim(dir))
! 		  if (.not. www) then
		  print *, "actual directory is:", trim(dir)
! 		  end if
		  dir_name=dir
		  skip = 1
	    case("--dir-global")
		  call getarg(i+1, dirglob%dir)
		  dirglob%valid = .true.
		  skip = 1            
	    case("")
		    CONTINUE
	    case default
		    print *, "incorrect command line options, your option was:"
		    print *, trim(arg)
		    ERROR STOP
	  end select
	else
	  do j=1, ubound(skip,1)
	    if (skip(j) == 1) then
	      skip(j) = 0
	      EXIT  
	    end if
	  end do
	  if (ubound(skip,1) > 1 .and. maxval(skip) == 0) then
	    deallocate(skip)
	    allocate(skip(1))
	  end if
	end if
      end do
      
      if (www) then
	call find_unit(fileid,1000)
	call system("rm -rf 4www")
	call system("mkdir 4www")
	open(unit=fileid, file="4www/mypid", action="write", status="replace")
	write(unit=fileid, fmt=*) getpid()
	close(fileid)
      end if
    
    end subroutine get_cmd_options

    subroutine set_global_vars()
      use globals
      use core_tools
      use pde_objs

      integer :: i_err
      
      
      select case (adjustl(trim(drutes_config%name)))
	case("RE_std", "RE_rot", "REstdH", "RErotH", "boussi", "ADEstd")
	  pde_common%processes = 1
	case("ADE_RE_std", "ADE_REstdH", "ADE_RE_rot", "ADE_RErotH", "ADEstd_kinsorb")
	  pde_common%processes = 2
	case("RE_mod", "REtest", "ADE_RE_std_kinsorb", "ADE_REstdH_kinsorb", "ADE_RE_rot_kinsorb", &
	     "ADE_RErotH_kinsorb")
	  pde_common%processes = 3
      end select

      backup_runs = 0

      !set the terminal output unit ID
      if (.not. (terminal_assigned)) then
	select case(print_level)
	  case(0)
	    terminal = 6
	  case(1,-1)
	    call find_unit(terminal, 2000)
	    if (print_level == 1) then
	      open(unit=terminal, file="out/screen.log", action="write", status="replace", iostat=i_err)
	      if (i_err /=0) then
		print *, "ERROR: check if directory out/ exist!"
		ERROR STOP
	      end if
	    else
	      open(unit=terminal, file="/dev/null", action="write", status="replace", iostat=i_err)
	      if (i_err /=0) then
		print *, "ERROR: file /dev/null does not exist, are you running GNU/Linux based os?"
		print *, "       if not do not set print level -1 in drutes.conf/global.conf"
		ERROR STOP
	      end if
	    end if
	end select
      end if

      solver_call = 0

      postpro_run = -1_ikind
      postpro_dec = 1_ikind
      

    end subroutine set_global_vars


    subroutine init_observe()
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use core_tools


      integer(kind=ikind) :: i, j, k, point, dec
      real(kind=rkind), dimension(:,:), allocatable :: domain
      logical :: error_stop = .false.

      allocate(domain(ubound(elements%data,2), drutes_config%dimen))

      observation_array(:)%element = 0
      
      do i=1, ubound(observation_array,1)
	allocate(observation_array(i)%cumflux(ubound(pde,1)))
	observation_array(i)%cumflux = 0.0_rkind
      end do
      

      do i = 1, ubound(observation_array,1)
	do j=1, elements%kolik
	  do k=1, ubound(elements%data,2)
	    domain(k,:) = nodes%data(elements%data(j,k),:)
	  end do
	  if (inside(domain, observation_array(i)%xyz)) then

	    observation_array(i)%element = j
	    EXIT
	  else
	    CONTINUE
	  end if
	end do
      end do


      do i=1, ubound(observation_array,1)
	if (observation_array(i)%element == 0) then
	  print *, "ERROR: observation point: ", i, "with coordinates: ", observation_array(i)%xyz, "lies out of domain"
	  error_stop = .true.
	end if
      end do

      if (error_stop) then
	ERROR STOP "error in observation points definition, check drutes.conf/global.conf"
      end if



      do point=1, ubound(observation_array,1)
	do i=1, ubound(pde,1)
	  call find_unit(pde(i)%obspt_unit(point),5000)
	  dec = 0
	  do 
	    dec = dec+1
	    if (point/(10**dec) < 1) then
	      EXIT
	    end if
	  end do
	  call open_obs_unit(pde(i), point, dec)
	end do
      end do

      deallocate(domain)
      

    end subroutine init_observe
    
    
    subroutine init_measured()
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use core_tools 
      
      integer(kind=ikind) :: i,j
      real(kind=rkind) :: act_dist
      
      
      act_dist = huge(act_dist)
      do i=1, ubound(measured_pts,1)
	do j=1, nodes%kolik
	  if (dist(measured_pts(i)%xyz, nodes%data(j,:)) < act_dist) then
	    measured_pts(i)%node = j
	    act_dist = dist(measured_pts(i)%xyz, nodes%data(j,:))
	    if (act_dist <= epsilon(act_dist)) EXIT
	  end if
	end do
      end do
      
      outer_boundaries = maxval(nodes%edge)
      
      do i=1, ubound(measured_pts,1)
       nodes%edge(measured_pts(i)%node) = outer_boundaries + i
      end do
      
	  

      
      
    end subroutine init_measured


    subroutine open_obs_unit(pde_loc, name, decimals)
      use typy
      use globals
      use global_objs
      use pde_objs

      class(pde_str), intent(in) :: pde_loc
      integer(kind=ikind), intent(in) :: name
      integer(kind=ikind), intent(in) :: decimals
      character(len=256) :: fileid
      character(len=64) :: forma
      character(len=10), dimension(3) :: xyz


      xyz(1) = "x"
      xyz(2) = "z"
      xyz(3) = "y"
      
      write(unit=forma, fmt="(a, I16, a)") "(a, a, a, I", decimals, ", a)"
      write(unit=fileid, fmt=forma) "out/obspt_", adjustl(trim(pde_loc%problem_name(1))), "-", name, ".out"
     
      if (.not. drutes_config%run_from_backup) then
	open(unit=pde_loc%obspt_unit(name), file=fileid, action="write", status="replace")
	
	write(unit=pde_loc%obspt_unit(name), fmt=*) "#        time                      ", &
	  trim(pde_loc%solution_name(2)), "            ", &
	!  trim(pde_loc%solution_name(2)), " (avg over el.)", &
	!  "       ", trim(pde_loc%mass_name(2)), "       ",&
	  trim(pde_loc%flux_name(2)), "   in    ", xyz(1:drutes_config%dimen), "     directions", "   cumulative flux"
	write(unit=pde_loc%obspt_unit(name), fmt=*) &
	  "#-----------------------------------------------------------------------------------------------"
	write(unit=pde_loc%obspt_unit(name), fmt=*)
      else
	open(unit=pde_loc%obspt_unit(name), file=fileid, action="write", access="append", status="old")
      end if

	  
    end subroutine open_obs_unit
    
    
    subroutine init_obstimes()
      use typy
      use globals
      use global_objs
      
      real(kind=rkind), dimension(:,:), allocatable :: tables
      real(kind=rkind) :: curtime, dt
      integer(kind=ikind) :: i, j, insert
      type(smartarray_int) :: loc
      
      
      allocate(tables(ubound(observe_time,1),2))      
      tables(:,1) = observe_time(1:ubound(tables,1))%value
      
      if (observe_info%anime) then
	dt = end_time/observe_info%nframes
	do i=1, observe_info%nframes
	  tables(i + ubound(observe_time,1) - observe_info%nframes,1) = i*dt
	end do	
      end if
 
      curtime = 0.0

      
      i=0
      do
	where (tables(:,1) - curtime > epsilon(curtime))
	  tables(:,2) = tables(:,1) - curtime
	else where
	  tables(:,2) = huge(dt)
	end where
	
	call loc%nrfill(1_ikind*minloc((tables(:,2)),1))
	
	do j=1, ubound(observe_time,1)
	  if (abs(tables(j,1) - tables(loc%data(1),1)) < epsilon(dt) .and. j/=loc%data(1)) then
	    call loc%nrfill(j)
	  end if
	end do
	
	do j=1, loc%pos
	  insert = loc%data(j)
	  i=min(i+1, ubound(observe_time,1))
	  
	  if (i == 0) then
	    RETURN
	  end if
	  
	  if (insert >  ubound(observe_time,1) - observe_info%nframes .and. observe_info%anime) then
	    observe_time(i)%name = "avi"
	  else
	    observe_time(i)%name = "obs"
	  end if
	  
	  observe_time(i)%value = tables(insert,1)
	  curtime = observe_time(i)%value
	end do

	call loc%clear()
	
	if (i==ubound(observe_time,1)) then
	  EXIT
	end if
	
      end do
      
    end subroutine init_obstimes


end module drutes_init