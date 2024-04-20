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

!> \file objfnc.f90
!! \brief Use for inverse modeling. 
!<

!>
!! Tool for inverse modeling.  User can specify file with experimental data, DRUtES can compute value of objective function
!<



module objfnc
  use typy
  
  private :: reader
  private :: get_objval_standard
  private :: read_model
  private :: time_exceeded
  
  !> object with a method definition for objective function computation
  type, public :: objval_obj
    !> when performing inverse analyses, by setting this value .true., drutes can compute value(s) of objective function(s)
    logical :: compute
    !> by setting this value if some defined code runtime exceeded, then drutes can take different actions according to CPUlim_method
    logical :: limit_CPU=.false.
    !> if CPUlim_method == E DRUtES exits computation, objective value is computed
    !! if CPUlim_method == P Picard criterion is updated to the value Picard_new
    !<
    character(len=1) :: CPUlim_method
    !> updated value of the Picard criterion, typically some huge number, so that the Picard process degenerates into semi-explicit method
    real(kind=rkind) :: Picard_new
    !> time in seconds of the limitting CPU time
    real(kind=rkind) :: CPUtime_max
    contains 
      procedure, nopass :: read_config=>reader
      procedure, nopass :: getval=>get_objval_standard   
      procedure, nopass :: toolong=>time_exceeded
  end type objval_obj  
  
  type(objval_obj), public :: objval
  
  integer, private, save :: exp_file
  
  type, private :: point_str
    real(kind=rkind), dimension(:), allocatable :: time
    real(kind=rkind), dimension(:,:), allocatable :: data
  end type point_str
  
  type, private :: ram_limit_str
    logical :: set
    real(kind=rkind) :: memsize=0.0
    character(len=2) :: units
  end type ram_limit_str

  type(point_str), dimension(:), allocatable, private, save :: exp_data
  type(point_str), dimension(:), allocatable, private, save :: model_data
  integer(kind=ikind), dimension(:), allocatable, private, save :: obs_ids
  integer(kind=ikind), dimension(:), allocatable, private, save :: noprop, pde_comp
  integer(kind=ikind), dimension(:,:), allocatable, private, save :: columns
  character(len=4096), private, save :: fileinputs
  type(ram_limit_str), private, save :: ram_limit
  integer(kind=ikind), private, save ::  no_pdes
  integer, dimension(:), allocatable, private, save :: datafiles
  
  integer, private, save :: fileid, expfile, ierr
  integer(kind=ikind), private :: n, i, counter, tmpbound, expcols, i1, j, low, top, skipcount, datacount, l
  character(len=4096), private :: msg
  character(len=4), private :: units
  real(kind=rkind), private :: r
  real(kind=rkind), dimension(:), allocatable, private :: tmpdata
  real(kind=rkind), private, save :: memsize, corr_memsize
  logical, dimension(:), allocatable, private, save :: skipid
  logical, private :: go4skip, processed
  
  integer(kind=ikind), private :: pos, k
      
  type, private :: errors_str
    real(kind=rkind), dimension(:), allocatable :: val
  end type
      
  type(errors_str), dimension(:), allocatable, private, save :: errors
  logical, private :: inlast 
  integer, private :: outfile
  real(kind=rkind), private :: dt, slope, modval, suma
  
  
  
  
  contains
    
    !>action to be taken if the CPU run time is longer than the allowed maximal time
    !! this function is linked to objval%toolong, replace with something more elegant if you like.
    !<
    subroutine time_exceeded(cpu_time)
      use typy
      use globals
      use core_tools
      
      !> current CPU time, must be real(4) because of old etime function definitions
      real(4) :: cpu_time
      logical :: written=.false.
      
      select case(objval%CPUlim_method)
        case("E")
          call write_log("Maximal CPU time for evaluating your objective function has exceeded, &
           DRUtES will stop, your objective function will be still evaluated, but some part of your &
            simulation will be missing")
          cpu_time_limit=.true.
          cpu_max_time = cpu_time - epsilon(cpu_time)
        case("P")
          if (.not. written) then
            call write_log("Maximal CPU time for evaluating your objective function has exceeded, &
            DRUtES will now increase the Picard criterion from the previous value:", real1=iter_criterion, &
            text2="for the new value", real2=objval%Picard_new)
            written=.true.
          end if
          
          iter_criterion=objval%Picard_new
      end select
      
    end subroutine time_exceeded
      
    
    subroutine reader()
      use typy
      use readtools
      use core_tools
      use globals
      use debug_tools
      use pde_objs
      

      
      open(newunit=fileid, file="drutes.conf/inverse_modeling/objfnc.conf", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("unable to open file drutes.conf/inverse_modeling/objfnc.conf with objective function &
        configuration, exiting DRUtES....")
        ERROR STOP
      end if
      
      write(msg, *) "Set the total number of components you want to model", new_line("a"), &
                    "    e.g. your objective function consists of concetration and pressure head values",  new_line("a"), &
                    "    then this number is equal 2"
                    
      call fileread(no_pdes, fileid, errmsg=msg, ranges=(/1_ikind, 1_ikind*ubound(pde,1)/))
      
      call read_sep(fileid)
      
      write(msg, *) "The number of PDE component is defined as follows:" , new_line("a"), &
           "      Richards equation - always 1, no other option",  new_line("a"), &
           "      Dual permeability problem - 1 for matrix, 2 for fracture",  new_line("a"), &
           "      Advection-dispersion-reaction equation with advection defined in conf files - always 1 , no other option", &
           new_line("a"), &
           "      Advection-dispersion-reaction equation with advection computed - 1 for e.g. pressure head, &
            2 for concentration" , & 
           new_line("a"), &
           "      Advection-dispersion-reaction equation with advection computed and kinetic sorption - & 
           1 for e.g. pressure head, 2 for &
           concentration in liquid phase, 3 for concentration in solid phase", new_line("a"), new_line("a"), &
           "I M P O R T A N T !!! the number of lines with pde component ids has to be equal to the number of &
           components defined above!!"
           
      
      allocate(pde_comp(no_pdes))
      
      do i=1, no_pdes 
        call fileread(pde_comp(i), fileid, errmsg=msg, ranges=(/1_ikind, 1_ikind*ubound(pde,1)/))
      end do
           
      call read_sep(fileid)
      
      write(msg, *) "Check the number of your points for constructing objective function, it should be equal or lower than ", &
          "the number of observation points and at least 1.", new_line("a"), &
          "   Your number of observation points is: ",   ubound(observation_array,1)
      
      call fileread(n, fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      
      call read_sep(fileid)
      
      
      allocate(obs_ids(n*no_pdes))
      
      msg="Are the numbers of observation points for evaluating your objective function correct?"
      
      do i=1, n
        call fileread(obs_ids(i), fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      end do

      call read_sep(fileid)
      
      allocate(noprop(n*no_pdes))
      
      write(msg, *) "   For each observation point you must specify number of properties you want to check, &
          the range is 1 till 4", &
        new_line("a"), "1st property is typically solution, 2nd is mass (e.g. water content), 3rd is flux, 4th cummulative flux.", & 
        new_line("a"), "See the head of the output file of the observation points!!" 
      do i=1,n        
        call fileread(noprop(i), fileid, ranges=(/1_ikind, 4_ikind/), errmsg=msg)
      end do 
      
      call read_sep(fileid)
      allocate(columns(ubound(noprop,1), (maxval(noprop))))
      
      write(msg, *) "Is the number of columns for evaluating your objective function correct?", new_line("a"), &
          " 1st column is reserved for time, start with 2nd column, which is typically reserved for the primary solution"
      
      columns = 0
      
      do i=1, ubound(noprop,1)
!       time, val, massval, advectval(1:D), observation_array(i)%cumflumemsizex(proc) - in total 4 properties + time
        call fileread(columns(i, 1:noprop(i)), fileid, ranges=(/2_ikind, 5_ikind/), errmsg=msg, checklen=.TRUE.)
      end do
      
      
      call read_sep(fileid)
      
      call fileread(fileinputs, fileid)  
      
      call read_sep(fileid)
      
      call write_log("The file with inputs for inverse modeling is: ", text2=cut(fileinputs))
      
      
      expcols=0
      
      do i=1, ubound(noprop,1)
        expcols=expcols+noprop(i)
      end do
      
            
      open(newunit=expfile, file=cut(fileinputs), status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "the file with your inputs doesn't exist"
        print *, "----------------"
        print *, "HINT: If this run is not a part of automatic calibration script,"
        print *, "and you want to avoid this stop, simply define [n] for Do inverse modeling"
        print *, "in file     d r u t e s . c o n f / g l o b a l . c o n f"
        print *, "----------------"
        print *, "----------------"
        ERROR STOP
      end if
      
      call fileread(ram_limit%set, fileid)
      
      call read_sep(fileid)
      
      call fileread(ram_limit%memsize, fileid)
      
      call read_sep(fileid)
      
      call fileread(ram_limit%units, fileid, options=(/"kB", "MB", "GB"/))
      
      call read_sep(fileid)
      
      if (ierr /= 0) then
        write(msg, *) "ERROR!, You have specified bad path for input file with experimental data for inverse modeling.", & 
        new_line("a"), &
        " The path you have specified is: ", cut(fileinputs), new_line("a"), &
        " The path should start with your DRUtES root directory", new_line("a"), &
        " e.g. drutes.conf/inverse_modeling/inputs.dat"
        call file_error(fileid, msg)
      end if
      
      call fileread(objval%limit_CPU, fileid)
      
      call read_sep(fileid)
      
      if (objval%limit_CPU) then
        
        call comment(fileid)

        read(unit=fileid, fmt=*, iostat=ierr) objval%CPUtime_max , units, objval%CPUlim_method
        
        if (ierr /= 0) then
          write(msg, fmt=*) "You have incorrect setup for limitting the CPU run time", new_line("a"), &
            "Set the following values:", new_line("a"), &
            "max. run time  |  unit (s , min, hrs, day) | action (E=exit the code, P=update & 
              Picard criterion)", new_line("a"), "See this example:", new_line("a"), &
              "100             min                      P"
          call file_error(fileid, errmsg=msg)
        end if
             
        
        select case(cut(units))
          case("s")
            CONTINUE
          case("min")
            objval%CPUtime_max = objval%CPUtime_max*60
          case("hrs")
            objval%CPUtime_max = objval%CPUtime_max*3600
          case("day")
            objval%CPUtime_max = objval%CPUtime_max*86400
          case default
            write(msg, *) "You have defined unsupported units for CPU time limit. the supported units are: s, min, hrs, day", &
              new_line("a"), "You have specified: ", cut(units), " , correct it!!"
            call file_error(fileid, errmsg=msg)
        end select
        
        if (objval%CPUlim_method=="P") then
          msg="The new update for the Picard criterion defined in drutes.conf/inverse_modeling/objfnc.conf &
                must be greater than the current Picard criterion defined in drutes.conf/global.conf &
                otherwise it is useless"
                
          call fileread(objval%Picard_new, fileid, &
          ranges=(/iter_criterion+epsilon(iter_criterion), huge(iter_criterion)/), errmsg=msg)
        end if
        
      end if 
      
      
      
      close(fileid)
      
      counter=0
      do 
        counter=counter+1
        call comment(expfile)
        read(unit=expfile, fmt=*, iostat=ierr) r
        if (ierr /= 0) then
          counter = counter-1
          EXIT
        end if
      end do
      
      

      allocate(exp_data(ubound(noprop,1)))
      
      allocate(exp_data(1)%time(counter))
      
      do i=1, ubound(noprop,1)
        allocate(exp_data(i)%data(counter, noprop(i)))
      end do
      
      
      close(expfile)
       
      open(newunit=expfile, file=cut(fileinputs), status="old", action="read")
      
      write(msg, *) "ERROR! You have either wrong definition in drutes.conf/inverse_modeling/objfnc.conf ", new_line("a") , &
        "or wrong number of columns in", trim(fileinputs), new_line("a"), &
        "e.g. for 2 observation points and 2 properties (in total 4 values per time) you need", new_line("a"), &
        "5 columns -> 1st col. = time, col. 2-5  = your properties"
        
        
      allocate(tmpdata(expcols+1))
      
      do i=1, counter-1     
        call fileread(tmpdata, expfile, checklen=.TRUE.)
        exp_data(1)%time(i) = tmpdata(1)
        low=2
        do j=1, ubound(exp_data,1)
          top=low+noprop(j)-1
          exp_data(j)%data(i, 1:noprop(j)) = tmpdata(low:top)
          low=top+1
        end do  
      end do
  

      deallocate(tmpdata)
      

    end subroutine reader
    
    
    subroutine read_model()
      use typy
      use readtools
      use core_tools
      use globals
      use debug_tools
      use pde_objs
    
      integer(kind=ikind) :: massdim
    
      allocate(model_data(ubound(obs_ids,1)*no_pdes))
      
      allocate(datafiles(ubound(obs_ids,1)*no_pdes))
      
      
      do i=1, no_pdes
        do j=1, ubound(obs_ids,1)
          open(newunit=datafiles((i-1)*ubound(obs_ids,1)+j), file=pde(pde_comp(i))%obspt_filename(obs_ids(j)), action="read", &
                 status="old", iostat=ierr)
          if (ierr/=0) then
            print *, "error opening files with observation points"
            print *, "this is a bug"
            print *, "called from objfnc::reader"
            print *, "contact Michal -> michalkuraz@gmail.com"
            ERROR STOP
          end if
        end do
      end do
      
      counter = 0
      
      do 
        counter=counter+1
        call comment(datafiles(1))
        read(unit=datafiles(1), fmt=*, iostat=ierr) r
        if (ierr /= 0) then
          counter = counter-1
          EXIT
        end if     
      end do
      
      memsize = rkind*(ubound(datafiles,1)+1)*counter
      
      select case(ram_limit%units)
        case("kB")
          ram_limit%memsize = ram_limit%memsize*1e3
        
        case("MB")
          ram_limit%memsize = ram_limit%memsize*1e6
        
        case("GB")
          ram_limit%memsize = ram_limit%memsize*1e9
        
      end select
      
      if (ram_limit%set .and. memsize > ram_limit%memsize) then
      
        corr_memsize = int(ram_limit%memsize/(rkind*no_pdes*ubound(obs_ids,1)))*(rkind*no_pdes*ubound(obs_ids,1))
        
        write(msg, *) "of RAM for objective function computation. The datafiles exceeded your RAM limit, and thus", &
                    100-int(corr_memsize/memsize*100),"% of your output data will be skipped."
        go4skip=.true.
        
      else
        corr_memsize = memsize
        write(msg, *) "of RAM for objective function computation."
        go4skip=.false.
      end if
      
      
      if (corr_memsize < 999) then
        call write_log("DRUtES will allocate: ", int1=1_ikind*nint(corr_memsize), text2="B of RAM", text3=trim(msg))
        else if (corr_memsize > 999 .and. corr_memsize < 1e6) then
          call write_log("DRUtES will allocate: ", real1=corr_memsize/1e3, text2="kB", text3=trim(msg))
        else if (corr_memsize > 1e6 .and. corr_memsize < 1e9 ) then
          call write_log("DRUtES will allocate: ", real1=corr_memsize/1e6, text2="MB", text3=trim(msg))
        else 
          call write_log("DRUtES will allocate: ", real1=corr_memsize/1e9, text2="GB", text3=trim(msg))
      end if
      
      close(datafiles(1))
      
      open(newunit=datafiles(1), file=pde(pde_comp(1))%obspt_filename(obs_ids(1)), action="read", status="old")
      
      if (go4skip) then
        
        allocate(skipid(counter))
        
        skipid = .false.
        
        skipcount =  (int((memsize-corr_memsize)/memsize*counter)+1)
        
        i=0
        do 
          call random_seed()
          n = int(counter*rand(0))
          if (n>0) then
            if (.not. skipid(n)) then
              i = i + 1
              skipid(n) = .true.
            end if
          end if
            
          if (i == skipcount) EXIT
        end do

        datacount = counter-skipcount
      else
        datacount = counter
      end if
      
      allocate(model_data(1)%time(datacount))
      
      do i=1, ubound(model_data,1)
        allocate(model_data(i)%data(datacount, noprop(i)))
      end do
      
    

      do i=1, ubound(model_data,1)
        if (allocated(tmpdata)) deallocate(tmpdata)
        
        if(pde(pde_comp(1))%print_mass) then
          massdim = ubound(pde(pde_comp(1))%mass,1)
        else
          massdim = 0
        end if
        
        allocate(tmpdata(3 + drutes_config%dimen + massdim))
        n=0
        do l=1, counter
          processed=.false.
          call fileread(tmpdata, datafiles(i))

          if (allocated(skipid)) then
            if (.not. skipid(l)) then
              processed = .true.
              n=n+1
            end if
          else
            processed = .true.
            n=n+1
          end if
          
          if (i==1 .and. processed) then
            model_data(i)%time(n) = tmpdata(1)
          end if
        
          if (processed) then
            do j=1, noprop(i)
              model_data(i)%data(n,j) = tmpdata(columns(i,j))
            end do
          end if
        end do
      end do  

    
    end subroutine read_model
    
    subroutine get_objval_standard()
      use typy
      use debug_tools
      
      
      call read_model()


      allocate(errors(ubound(model_data,1)))
      
      do i=1, ubound(errors,1)
        allocate(errors(i)%val(noprop(i)))
        errors(i)%val = 0
      end do
      
      
      pos = 1
      do j=1, ubound(exp_data(1)%time,1) 
        do k=pos, ubound(model_data(1)%time,1)-1
          if (exp_data(1)%time(j) >=  model_data(1)%time(k) .and. exp_data(1)%time(j) < model_data(1)%time(k+1)) then
            dt = model_data(1)%time(k+1) - exp_data(1)%time(j)
            if (dt < model_data(1)%time(k+1)*epsilon(dt)) then
              inlast = .true.
            else
              inlast = .false.
            end if
            
            pos = k
            
            n=0
            do i=1, ubound(errors,1)
              do l=1, ubound(errors(i)%val,1)
                n=n+1
                if (.not. inlast) then
                  slope = (model_data(n)%data(pos+1,l) - model_data(n)%data(pos,l))/ &
                          (model_data(1)%time(pos+1) - model_data(1)%time(pos))
                  modval = model_data(n)%data(pos+1,l) - slope*dt
                else
                  modval = model_data(n)%data(pos+1,l)
                end if
                
                errors(i)%val(l) = errors(i)%val(l) + (modval - exp_data(n)%data(j,l))*(modval - exp_data(n)%data(j,l))
                
              end do
            end do
           
            EXIT
          end if
        end do
      
      end do
        
        
   
      do i=1, ubound(errors,1)  
        errors(i)%val = sqrt(errors(i)%val/pos)
      end do
      
      open(newunit=outfile, file="out/objfnc.val", status="new", action="write")
      
      write(outfile, *) "# values of objective functions"
      
      do i=1, ubound(errors,1)
        do j=1, ubound(errors(i)%val,1)
          write(outfile, *) errors(i)%val(j)
        end do
      end do
      
      
      close(outfile)
                                
    end subroutine get_objval_standard
  


end module objfnc
