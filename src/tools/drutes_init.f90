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

!> \file drutes_init.f90
!! \brief DRUtES init procedures.
!<

module drutes_init
  public  :: parse_globals, init_observe, init_measured, get_cmd_options, pde_constructor, init_bcfluxes
  private :: set_global_vars,  open_obs_unit, init_obstimes


  contains
    !> 2nd procedure from main
    !! Procedure reading drutes.conf/global.conf file, opens out/DRUtES.log, out/cgit.count
    !<
    subroutine parse_globals()
      use typy
      use globals
      use pde_objs
      use core_tools
      use read_inputs
      use readtools
      use debug_tools
      use printtools
      use global4solver


      integer :: i_err, i
      character(len=256) ::  writer
      character(len=8192) :: dirname
      real(kind=rkind) :: const


      open(newunit=logfile, file="out/DRUtES.log", action="write", status="replace", iostat=i_err)
      
      if (i_err /= 0) then
        i_err=system("mkdir out")
        open(unit=logfile, file="out/DRUtES.log", action="write", status="replace", iostat=i_err)
        if (i_err /= 0) then
          print *, "unable to open and create directory out (is it UN*X type os?), called  from drutes_init::parse_globals"
          ERROR STOP
        end if
	  
      end if
      
      call print_logo(logfile)
    
      write(unit=logfile, fmt=*, iostat = i_err) "DRUtES is initializing, reading input files"

      open(newunit=file_itcg, file="out/itsolver.count", action="write", status="replace", iostat=i_err)
      
      open(newunit=file_picard, file="out/picard_it.count", action="write", status="replace", iostat=i_err)
      
      call write_log(text="total number of parallel images is:", int1=1_ikind*NUM_IMAGES())

      call write_log("reading input files and allocating...")

      open(newunit=file_global,file="drutes.conf/global.conf", action="read", status="old", iostat = i_err)

      if (i_err /= 0) then
        print *, "missing drutes.conf/global.conf file"
        ERROR STOP
      end if
      
      open(newunit=file_solver,file="drutes.conf/solver.conf", action="read", status="old", iostat = i_err)
      
      
      if (i_err /= 0) then
        print *, "missing drutes.conf/solver.conf file"
        ERROR STOP
      end if
      
      call read_global()
      
      call read_solverconfig()
      
      if (record_solver_time) open(newunit=solver_time_file,file="out/solver.time", action="write", status="replace")
     
      call set_global_vars()
      
      call init_obstimes()

       !read mesh data

      select case(drutes_config%dimen)
        case(1)
          open(newunit=file_mesh, file="drutes.conf/mesh/drumesh1d.conf", action="read", status="old", iostat=i_err)
          if (i_err /= 0) then
            print *, "missing drutes.conf/mesh/drumesh1d.conf file"
            ERROR STOP
          end if
          call read_1dmesh_int()

        case(2)
          select case(drutes_config%mesh_type)
            case(1)
              open(newunit=file_mesh, file="drutes.conf/mesh/drumesh2d.conf", action="read", status="old", iostat=i_err)
              if (i_err /= 0) then
                print *, "missing drutes.conf/mesh/drumesh2d.conf file"
                ERROR STOP
              end if
              call read_2dmesh_int()
            case(2)
              open(newunit=file_mesh, file="drutes.conf/mesh/mesh.t3d", action="read", status="old", iostat=i_err)
              if (i_err /= 0) then
                print *, "missing drutes.conf/mesh/mesh.t3d file"
                ERROR STOP
              end if
              call read_2dmesh_t3d()
            case(3)
              open(newunit=file_mesh, file="drutes.conf/mesh/mesh.msh", action="read", status="old", iostat=i_err)
              if (i_err /= 0) then
                print *, "missing drutes.conf/mesh/mesh.msh file"
                ERROR STOP
              end if
              call read_2dmesh_gmsh()
            case(4)
              call read_ArcGIS()
            end select
        
          case(3)
              open(newunit=file_mesh, file="drutes.conf/mesh/mesh3D.msh", action="read", status="old", iostat=i_err)
              if (i_err /= 0) then
                print *, "missing drutes.conf/mesh/mesh3D.msh file"
                ERROR STOP
              end if
              call read_3dmesh_gmsh()         
            
        case default
          write(unit=terminal, fmt=*)"ERROR: unsupported problem dimension, the specified dimension was: ", drutes_config%dimen
          write(unit=terminal, fmt=*)"currently only 1D and 2D is supported"
          ERROR STOP 
      end select
      

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
      allocate(pde(i)%mass(1))
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
        pde(i)%mass(1)%val => dummy_scalar
        pde(i)%flux => dummy_vector
        pde(i)%dt_check => time_check_ok
        pde(i)%process_change => do_nothing
        pde(i)%getval => getvalp1
        pde(i)%order = i
      end do


    end subroutine pde_constructor
    
    !> 1st procedure in main
    !! analyzes options passed to drutes executable
    !! the following options are available
    !! -o optim : batch executions of drutes, in case of divergence of Picard method the Picard criterion is increased to "huge" value, so that it becomes semi explicit method
    !! --print-level 0 : default option everything prints to terminal, so output is not saved
    !! --print-level 1 : all terminal ouputs are saved in out/screen.logfile
    !! --print-level 2 : all terminal ouputs go to /dev/null 
    !! -v : only prints drutes version
    !! --tmax : user can specify maximal CPU time for code execution e.g. --tmax 5 min , available units are: s , min , hrs , day 
    !! --dir-local and --dir-global and www was supposed to be used for online preprocessor and postprocessor, but it hasn't been fully implemented yet, ignore it for now
    !<
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
            case("-h")
              print *, "Available switches are:"
              print *, "-o optim : if nonlinear solver fails to converge within minimal time step, the setup is changed into"
              print *, "           semiexplicit method. It means the nonlinearities are computed based on previous time step values"
              print *, "--print-level 1 : all screen output goes to out/screen.log"
              print *, "--print-level 2 : all screen output goes to /dev/null (use on Linux only)"
              print *, "-v prints DRUtES version"
              print *, "--tmax [value] [unit] : will limit CPU runtime to a given value. Available units are [s, min, hrs, day]"
              print *, "                        example:"
              print *, "                                bin/drutes --tmax 30 min"
              print *, "                        if the computation will not finish in 30 mins the time integration loop will exit"
              print *, "                        before reaching the end time"
              STOP
            case("-o")
              call getarg(i+1,oarg)
              select case(trim(oarg))
                case("www")
                  www=.true.
                case("optim")
                  optim=.true.
                case default
                  print *, "unrecognized option after -o , exiting...."
                  print *, "available options are:"
                  print *, " -o optim"
                  print *, "try bin/drutes -h"
                  ERROR STOP
              end select

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
                  print *, "try bin/drutes -h"
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
    
    
    end subroutine get_cmd_options

    subroutine set_global_vars()
      use globals
      use core_tools
      use pde_objs

      integer :: i_err
      

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
      use debug_tools

      integer(kind=ikind) :: i, j, k, point, dec
      real(kind=rkind), dimension(:,:), allocatable :: domain
      real(kind=rkind), dimension(:,:,:), allocatable :: domain3D
      logical :: error_stop = .false., foundel
      integer(kind=ikind) :: ord1, ord2, ord3, ii

      if (drutes_config%dimen < 3) then 
        allocate(domain(ubound(elements%data,2), drutes_config%dimen))
      else
        allocate(domain3D(4, 3, 3))
      end if

      observation_array(:)%element = 0
      
      do i=1, ubound(observation_array,1)
        allocate(observation_array(i)%cumflux(ubound(pde,1)))
        observation_array(i)%cumflux = 0.0_rkind
      end do
      

      do i = 1, ubound(observation_array,1)
        do j=1, elements%kolik

          if (drutes_config%dimen < 3) then
            do k=1, ubound(elements%data,2)
              domain(k,:) = nodes%data(elements%data(j,k),:)
            end do
            foundel = inside(domain, observation_array(i)%xyz)
          else
            do k=1, ubound(elements%data,2)
              select case(k)
                case(1)
                  ord1 = 1
                  ord2 = 2
                  ord3 = 3
                 case(2)
                  ord1 = 2
                  ord2 = 3
                  ord3 = 4
                 case(3)
                  ord1 = 1
                  ord2 = 3
                  ord3 = 4
                 case(4)
                  ord1 = 1
                  ord2 = 2
                  ord3 = 4
                end select
              
                domain3D(k, 1, :) = nodes%data(elements%data(j, ord1),:)
                domain3D(k, 2, :) = nodes%data(elements%data(j, ord2),:)
                domain3D(k, 3, :) = nodes%data(elements%data(j, ord3),:)
              end do
              
              foundel = inside3D(domain3D, observation_array(i)%xyz)

            end if
          
            
          if (foundel) then
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

      
      do i=1, ubound(pde,1)
        allocate(pde(i)%obspt_filename(ubound(observation_array,1)))
      end do


      do point=1, ubound(observation_array,1)
        do i=1, ubound(pde,1)
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

    end subroutine init_observe
    
    
    subroutine init_bcfluxes()
      use typy
      use globals
      use global_objs
      use debug_tools
      use geom_tools
      
      integer(kind=ikind), dimension(:), allocatable :: edgeid
      integer(kind=ikind), dimension(:), allocatable :: nodeid
      integer(kind=ikind) :: ndwrt, el, i, j, nd, D, edge2wrt, extnd, ii, bc
      integer(kind=ikind), dimension(:), allocatable :: sizes, elcounter
      character(len=12) :: filename
      logical :: found
      
      D = drutes_config%dimen
            
      
      allocate(bcfluxes(101:maxval(nodes%edge)))
      allocate(sizes(101:maxval(nodes%edge)))
      allocate(elcounter(101:maxval(nodes%edge))) 
      
      allocate(edgeid(ubound(elements%data,2)))
      allocate(nodeid(D))
      
      sizes = 0
      do i=1, nodes%kolik
        if (nodes%edge(i) >= 101) then
          sizes(nodes%edge(i)) = sizes(nodes%edge(i)) + 1
        end if
      end do
	  
      sizes = sizes * D

      do i=101, ubound(bcfluxes,1)
        allocate(bcfluxes(i)%bcel(sizes(i)))
      end do

      elcounter = 0
      do i=1, elements%kolik
        edgeid = 0
        ndwrt = 0
        do j=1, ubound(elements%data,2)
          edgeid(j) = nodes%edge(elements%data(i,j))
          if (edgeid(j) > 100) then
            ndwrt = ndwrt + 1
          end if
        end do
        
        found = .false.
        if (ndwrt >= D) then
          select case(D)
            case(3)
                ! 1 2 3
                ! 1 3 4
                ! 2 3 4
                ! 1 2 4
              if (edgeid(1) == edgeid(2) .and. edgeid(1) == edgeid(3)) then
                found=.true.
                edge2wrt = edgeid(1)
                extnd = elements%data(i,4)
                nodeid(:) =  elements%data(i,1:3)
              end if
            
              if (edgeid(1) == edgeid(3) .and. edgeid(1) == edgeid(4)) then
                found=.true.
                edge2wrt = edgeid(1)
                extnd = elements%data(i,2)
                nodeid(:) =  elements%data(i,[1,3,4])
              end if
              
              if (edgeid(2) == edgeid(3) .and. edgeid(2) == edgeid(4)) then
                found=.true.
                edge2wrt = edgeid(2)
                extnd = elements%data(i,1)
                nodeid(:) =  elements%data(i,2:4)
              end if
              
              if (edgeid(1) == edgeid(2) .and. edgeid(1) == edgeid(4)) then
                found=.true.
                edge2wrt = edgeid(1)
                extnd = elements%data(i,3)
                nodeid(:) =  elements%data(i,[1,2,4])
              end if
            case(2)
              ! 1 2
              ! 1 3
              ! 2 3
              if (edgeid(1) == edgeid(2)) then
                found = .true.
                edge2wrt = edgeid(1)
                extnd = elements%data(i,3)
                nodeid(:) =  elements%data(i,1:2)
              end if
              
              if (edgeid(1) == edgeid(3)) then
                found = .true.
                edge2wrt = edgeid(1)
                extnd = elements%data(i,2)
                nodeid(:) =  elements%data(i,[1,3])
              end if
              
              if (edgeid(2) == edgeid(3)) then 
                found = .true.
                edge2wrt = edgeid(2)
                extnd = elements%data(i,1)
                nodeid(:) =  elements%data(i,[2,3])
              end if
          end select
        end if
        
        if (found) then

          elcounter(edge2wrt) = elcounter(edge2wrt) + 1
          
          ii = elcounter(edge2wrt)
 
          bcfluxes(edge2wrt)%pos = ii
          bcfluxes(edge2wrt)%bcel(ii)%element = i
          bcfluxes(edge2wrt)%bcel(ii)%surfnode(1:D) = nodeid(1:D)
          bcfluxes(edge2wrt)%bcel(ii)%extnode = extnd
          
          select case(D)
            case(2)
              bcfluxes(edge2wrt)%bcel(ii)%n_out(1:D) = get_normals2D(nodes%data(nodeid(1),:), nodes%data(nodeid(2),:), &
                                                                     nodes%data(extnd,:))
              bcfluxes(edge2wrt)%bcel(ii)%area = dist(nodes%data(nodeid(1),:), nodes%data(nodeid(2),:))
              
            case(3)
              bcfluxes(edge2wrt)%bcel(ii)%n_out(1:D) = get_normals3D(nodes%data(nodeid(1),:), nodes%data(nodeid(2),:), &
                                                                      nodes%data(nodeid(3),:), nodes%data(extnd,:))
!              print *, bcfluxes(edge2wrt)%bcel(ii)%n_out(1:D)
              bcfluxes(edge2wrt)%bcel(ii)%area = triarea(nodes%data(nodeid(1),:), nodes%data(nodeid(2),:), nodes%data(nodeid(3),:))                                                     
          end select
        end if
      end do


      do bc=lbound(bcfluxes,1), ubound(bcfluxes,1)
        write(filename, fmt="(a,I3, a)") "out/", bc, ".flux"
        open(newunit=bcfluxes(bc)%fileid, file=filename, action="write", status="replace")
        write(bcfluxes(bc)%fileid, fmt=*) "# time,        flux,            cumulative flux"
        call flush(bcfluxes(bc)%fileid)
      end do


    end subroutine init_bcfluxes  
    
!    subroutine init_bcfluxes_old()
!      use typy
!      use globals
!      use global_objs
!      use debug_tools
      
!      integer(kind=ikind) :: i, j, nd1, nd2, elid, e, bc
!      character(len=12) :: filename
      
!      allocate(bcfluxes(101:maxval(nodes%edge)))
      
!      print *, "xasca" ; stop
      
!      do i=1, nodes%kolik
!        if (nodes%edge(i) > 0) then
!          do e=1, nodes%element(i)%pos
!            elid = nodes%element(i)%data(e)
!            nd1 = i
!            do j=1, ubound(elements%data,2)
!              if (elements%data(elid,j) /= nd1) then
!                nd2 = elements%data(elid,j)
!                if (nodes%edge(nd2) == nodes%edge(nd1)) then
!                  call bcfluxes(nodes%edge(nd1))%elements%nrfill(elid)
!                  call bcfluxes(nodes%edge(nd1))%nodes%nrfill(nd1)
!                  call bcfluxes(nodes%edge(nd1))%nodes%nrfill(nd2)
!                end if
!              end if
!            end do
!          end do
!        end if
!      end do
      
!      do bc=lbound(bcfluxes,1), ubound(bcfluxes,1)
!        write(filename, fmt="(a,I3, a)") "out/", bc, ".flux"
!        open(newunit=bcfluxes(bc)%fileid, file=filename, action="write", status="replace")
!        write(bcfluxes(bc)%fileid, fmt=*) "# time,        flux,            cumulative flux"
!      end do
      
    
!    end subroutine init_bcfluxes_old
    
    
    subroutine init_measured()
      use typy
      use globals
      use global_objs
      use pde_objs
      use geom_tools
      use core_tools 
      use debug_tools
      
      integer(kind=ikind) :: i,j
      real(kind=rkind) :: act_dist
      
      
      
      do i=1, ubound(measured_pts,1)
        act_dist = huge(act_dist)
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
      use printtools
      use debug_tools
      use core_tools

      class(pde_str), intent(in out) :: pde_loc
      integer(kind=ikind), intent(in) :: name
      integer(kind=ikind), intent(in) :: decimals
      character(len=64) :: forma
      character(len=10), dimension(3) :: xyz
      integer(kind=ikind) :: i, j
      character(len=512), dimension(:), allocatable :: fluxname


      xyz(1) = "x"
      xyz(2) = "z"
      
      if (drutes_config%dimen == 3) then
        xyz(2) = "y"
        xyz(3) = "z"
      end if
      
      write(unit=forma, fmt="(a, I16, a)") "(a, a, a, I", decimals, ", a)"
      write(unit=pde_loc%obspt_filename(name), fmt=forma) "out/obspt_", adjustl(trim(pde_loc%problem_name(1))), "-", name, ".out"
      
      if (allocated(pde_loc%fluxes)) then
        allocate(fluxname(ubound(pde_loc%fluxes,1)))
        do i=1, ubound(fluxname,1)
          write(fluxname(i), fmt=*) cut(pde_loc%fluxes(i)%name(2)),  " for: ", xyz(1:drutes_config%dimen), "     directions", &
            "   cumulative flux", "|----| "
        end do
      end if
      
     
      if (.not. drutes_config%run_from_backup) then
        if (ubound(pde_loc%mass_name,1) > 0) then
          open(newunit=pde_loc%obspt_unit(name), file=adjustl(trim(pde_loc%obspt_filename(name))), action="write", status="replace")
          
          if (cut(observe_info%fmt) == "pure") call print_logo(pde_loc%obspt_unit(name))
          
          
            i=1
            if (.not. allocated(pde_loc%fluxes)) then
              write(unit=pde_loc%obspt_unit(name), fmt=*) "#        time                      ", &
                trim(pde_loc%solution_name(2)), "            ", &
               "       ",  (/ (trim(pde_loc%mass_name(i,2)), i=1,ubound(pde_loc%mass_name,1) )   /), "       ",&
                trim(pde_loc%flux_name(2)), "   in    ", xyz(1:drutes_config%dimen), "     directions", "   cumulative flux"
              write(unit=pde_loc%obspt_unit(name), fmt=*) &
                "#-----------------------------------------------------------------------------------------------"
              write(unit=pde_loc%obspt_unit(name), fmt=*)
            else
              i = 1
              write(unit=pde_loc%obspt_unit(name), fmt=*) "#        time                      ", &
                trim(pde_loc%solution_name(2)), "            ", &
               "       ",  (/ (trim(pde_loc%mass_name(i,2)), i=1,ubound(pde_loc%mass_name,1) )   /), "       ",&
                (/ (cut(fluxname(i)),  i=1, ubound(fluxname,1) ) /)
              write(unit=pde_loc%obspt_unit(name), fmt=*) &
                "#-----------------------------------------------------------------------------------------------"
              write(unit=pde_loc%obspt_unit(name), fmt=*)
            end if
            
        else
          open(newunit=pde_loc%obspt_unit(name), file=adjustl(trim(pde_loc%obspt_filename(name))), action="write", status="replace")
          
          if (cut(observe_info%fmt) == "pure") call print_logo(pde_loc%obspt_unit(name))
          
          if (.not. allocated(pde_loc%fluxes)) then
            write(unit=pde_loc%obspt_unit(name), fmt=*) "#        time                      ", &
              trim(pde_loc%solution_name(2)), "            ", &
              trim(pde_loc%flux_name(2)), "   in    ", xyz(1:drutes_config%dimen), "     directions", "   cumulative flux"
            write(unit=pde_loc%obspt_unit(name), fmt=*) &
              "#-----------------------------------------------------------------------------------------------"
            write(unit=pde_loc%obspt_unit(name), fmt=*)
          else
            i=1
            write(unit=pde_loc%obspt_unit(name), fmt=*) "#        time                      ", &
              trim(pde_loc%solution_name(2)), "            ", &
             (/ (cut(fluxname(i)),  i=1, ubound(fluxname,1) ) /)
            write(unit=pde_loc%obspt_unit(name), fmt=*) &
              "#-----------------------------------------------------------------------------------------------"
            write(unit=pde_loc%obspt_unit(name), fmt=*)
          end if
        end if


      else
        open(newunit=pde_loc%obspt_unit(name), file=adjustl(trim(pde_loc%obspt_filename(name))), &
                    action="write", access="append", status="old")
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
