module objfnc
  use typy
  
  private :: reader
  public :: get_objval
  
  integer, private :: exp_file
  
  type, private :: point_str
    real(kind=rkind), dimension(:), allocatable :: time
    real(kind=rkind), dimension(:,:), allocatable :: data
  end type point_str
  
  type, private :: ram_limit_str
    logical :: set
    real(kind=rkind) :: memsize
    character(len=2) :: units
  end type ram_limit_str

  type(point_str), dimension(:), allocatable, private :: exp_data
  type(point_str), dimension(:), allocatable, private :: model_data
  integer(kind=ikind), dimension(:), allocatable, private :: obs_ids
  integer(kind=ikind), dimension(:), allocatable, private :: noprop
  integer(kind=ikind), dimension(:,:), allocatable, private :: columns
  character(len=4096), private :: fileinputs
  type(ram_limit_str), private :: ram_limit
  integer(kind=ikind), private :: pde_component, no_pdes
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
      
      integer(kind=ikind) :: n, i, counter, tmpbound, expcols, i1, j, low, top
      character(len=4096) :: msg
      real(kind=rkind) :: r
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      
      open(newunit=fileid, file="drutes.conf/inverse_modeling/objfnc.conf", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("unable to open file drutes.conf/inverse_modeling/objfnc.conf with objective function &
        configuration, exiting DRUtES....")
        ERROR STOP
      end if
      
      write(msg, *) "Set the total number of components you want to model correctly", new_line("a"), &
                    "    e.g. your objective function consists of concetration and pressure head values",  new_line("a"), &
                    "    then this numebr is equal 2"
                    
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
           concentration in liquid phase, 3 for concentration in solid phase"
           
      call fileread(pde_component, fileid, errmsg=msg, ranges=(/1_ikind, 1_ikind*ubound(pde,1)/))
           
      
      
      write(msg, *) "Check the number of your points for constructing objective function, it should be equal or lower than ", &
          "the number of observation points and at least 1.", new_line("a"), &
          "   Your number of observation points is: ",   ubound(observation_array,1)
      
      call fileread(n, fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      
      allocate(obs_ids(n))
      
      msg="Are the numbers of observation points for evaluating your objective function correct?"
      
      do i=1, n
        call fileread(obs_ids(i), fileid, ranges=(/1_ikind, 1_ikind*ubound(observation_array,1)/), errmsg=msg)
      end do

      allocate(noprop(n))
      
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
!       time, val, massval, advectval(1:D), observation_array(i)%cumflux(proc) - in total 4 properties + time
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
      
!       memsize = rkind*ubound(datafiles,1
!       
!       call write_log("DRUtES will allocate: ", int1=rkind*
      
    
    end subroutine reader
    
    subroutine get_objval()
      use typy
      
      call reader()
    
    end subroutine get_objval
  


end module objfnc
