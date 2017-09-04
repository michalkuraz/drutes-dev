module objfnc
  use typy
  
  private :: reader
  private :: get_objval_standard
  private :: read_model
  
  type, public :: objval_obj
    contains 
      procedure, nopass :: read_config=>reader
      procedure, nopass :: getval=>get_objval_standard   
  end type objval_obj  
  
  type(objval_obj), public :: objval
  
  integer, private :: exp_file
  
  type, private :: point_str
    real(kind=rkind), dimension(:), allocatable :: time
    real(kind=rkind), dimension(:,:), allocatable :: data
  end type point_str
  
  type, private :: ram_limit_str
    logical :: set
    real(kind=rkind) :: memsize=0.0
    character(len=2) :: units
  end type ram_limit_str

  type(point_str), dimension(:), allocatable, private :: exp_data
  type(point_str), dimension(:), allocatable, private :: model_data
  integer(kind=ikind), dimension(:), allocatable, private :: obs_ids
  integer(kind=ikind), dimension(:), allocatable, private :: noprop, pde_comp
  integer(kind=ikind), dimension(:,:), allocatable, private :: columns
  character(len=4096), private :: fileinputs
  type(ram_limit_str), private :: ram_limit
  integer(kind=ikind), private ::  no_pdes
  integer, dimension(:), allocatable :: datafiles
  
  
  
  
  contains
    
    subroutine reader()
      use typy
      use readtools
      use core_tools
      use globals
      use debug_tools
      use pde_objs
      
      integer :: fileid, expfile, ierr
      
      integer(kind=ikind) :: n, i, counter, tmpbound, expcols, i1, j, low, top, skipcount, datacount, l
      character(len=4096) :: msg
      real(kind=rkind) :: r
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      real(kind=rkind) :: memsize, corr_memsize
      logical, dimension(:), allocatable :: skipid
      logical :: go4skip, processed
      
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
           
      
      
      write(msg, *) "Check the number of your points for constructing objective function, it should be equal or lower than ", &
          "the number of observation points and at least 1.", new_line("a"), &
          "   Your number of observation points is: ",   ubound(observation_array,1)
      
      call fileread(n, fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      
      allocate(obs_ids(n*no_pdes))
      
      msg="Are the numbers of observation points for evaluating your objective function correct?"
      
      do i=1, n
        call fileread(obs_ids(i), fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      end do

      allocate(noprop(n*no_pdes))
      
      write(msg, *) "   For each observation point you must specify number of properties you want to check, &
          the range is 1 till 4", &
        new_line("a"), "1st property is typically solution, 2nd is mass (e.g. water content), 3rd is flux, 4th cummulative flux.", & 
        new_line("a"), "See the head of the output file of the observation points!!" 
      do i=1,n        
        call fileread(noprop(i), fileid, ranges=(/1_ikind, 4_ikind/), errmsg=msg)
      end do 
      
      
      allocate(columns(ubound(noprop,1), (maxval(noprop))))
      
      write(msg, *) "Is the number of columns for evaluating your objective function correct?", new_line("a"), &
          " 1st column is reserved for time, start with 2nd column, which is typically reserved for the primary solution"
      
      columns = 0
      
      do i=1, ubound(noprop,1)
!       time, val, massval, advectval(1:D), observation_array(i)%cumflumemsizex(proc) - in total 4 properties + time
        call fileread(columns(i, 1:noprop(i)), fileid, ranges=(/2_ikind, 5_ikind/), errmsg=msg, checklen=.TRUE.)
      end do
      
      
      call fileread(fileinputs, fileid)  
      
      call write_log("The file with inputs for inverse modeling is: ", text2=adjustl(trim(fileinputs)))
      
      expcols=0
      
      do i=1, ubound(noprop,1)
        expcols=expcols+noprop(i)
      end do
      
            
      open(newunit=expfile, file=adjustl(trim(fileinputs)), status="old", action="read", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "the file with your inputs doesn't exist"
        ERROR STOP
      end if
      
      call fileread(ram_limit%set, fileid)
      
      call fileread(ram_limit%memsize, fileid)
      
      call fileread(ram_limit%units, fileid, options=(/"kB", "MB", "GB"/))
      
      if (ierr /= 0) then
        write(msg, *) "ERROR!, You have specified bad path for input file with experimental data for inverse modeling.", & 
        new_line("a"), &
        " The path you have spefied is: ", adjustl(trim(fileinputs)), new_line("a"), &
        " The path should start with your DRUtES root directory", new_line("a"), &
        " e.g. drutes.conf/inverse_modeling/inputs.dat"
        call file_error(fileid, msg)
      end if
      
      
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
       
      open(newunit=expfile, file=adjustl(trim(fileinputs)))
      
      write(msg, *) "ERROR! You have either wrong definition in drutes.conf/inverse_modeling/objfnc.conf ", new_line("a") , &
        "or wrong number of columns in", trim(fileinputs), new_line("a"), &
        "e.g. for 2 observation points and 2 properties (in total 4 values per time) you need", new_line("a"), &
        "5 columns -> 1st col. = time, col. 2-5  = your properties"
        
        
      allocate(tmpdata(expcols+1))
      
      do i=1, counter     
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
      
      
      allocate(model_data(ubound(obs_ids,1)*no_pdes))
      
      allocate(datafiles(ubound(obs_ids,1)*no_pdes))
      
      
      do i=1, no_pdes
        do j=1, ubound(obs_ids,1)
          open(newunit=datafiles((i-1)*ubound(obs_ids,1)+j), file=pde(i)%obspt_filename(j), action="read", status="old", & 
               iostat=ierr)
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
      
      open(newunit=datafiles(1), file=pde(1)%obspt_filename(1), action="read", status="old")
      
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
      
      allocate(tmpdata(4+drutes_config%dimen))

      do i=1, ubound(model_data,1)
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
    
    end subroutine reader
    
    
    subroutine read_model()
    
    end subroutine read_model
    
    subroutine get_objval_standard()
      use typy
      use debug_tools
      
      integer(kind=ikind) :: i,j, pos, k, n, l
      
      type :: errors_str
        real(kind=rkind), dimension(:), allocatable :: val
      end type
      
      type(errors_str), dimension(:), allocatable :: errors
      logical :: inlast 
      integer :: outfile
      real(kind=rkind) :: dt, slope, modval, suma
      
      
      call reader()
      

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
        errors(i)%val = sqrt(errors(i)%val)
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
