module objfnc
  use typy
  
  private :: reader
  public :: get_objval
  
  integer, private :: exp_file
  
  type, private :: point_str
    real(kind=rkind), dimension(:), allocatable :: time
    real(kind=rkind), dimension(:), allocatable :: data
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
  
  
  
  contains
    
    subroutine reader()
      use typy
      use readtools
      use core_tools
      
      integer :: fileid, expfile, ierr
      
      integer(kind=ikind) :: n, i, counter, tmpbound, expcols, i1, i2
      character(len=4096) :: msg
      real(kind=rkind) :: r
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      
      open(newunit=fileid, file="drutes.conf/inverse_modeling/objfnc.conf", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        call write_log("unable to open file drutes.conf/inverse_modeling/objfnc.conf with objective function &
        configuration, exiting DRUtES....")
        ERROR STOP
      end if
      
      call fileread(n, fileid, ranges=(/0_ikind, huge(n)/))
      
      allocate(obs_ids(n))
      
      msg="Are the numbers of observation points for evaluating your objective function correct?"
      
      do i=1, n
        call fileread(obs_ids(i), fileid, ranges=(/1_ikind, ubound(observation_array,1)/), errmsg=msg)
      end do

      allocate(noprop(n))
      
      msg="For each observation point you must specify number of properties you want to check, the range is 1 till 4"
      do i=1,n        
        call fileread(noprop(i), fileid, ranges=(/1_ikind, 4_ikind/), errmsg=msg)
      end do  
      
      allocate(columns(ubound(noprop,1), (maxval(noprop))))
      
      msg="Is the number of properties for evaluating your objective function correct?"
      
      columns = 0
      
      do i=1, ubound(noprop,1)
!       time, val, massval, advectval(1:D), observation_array(i)%cumflux(proc) - in total 4 properties + time
        call fileread(columns(i, 1:noprop(i)), fileid, ranges=(/0_ikind, 4_ikind/), errmsg=msg, checklen=.TRUE.)
      end do
      
      call fileread(fileinputs, fileid)
      
      call write_log("The file with inputs for inverse modeling is", text2=trim(fileinputs))
      
      expcols=0
      
      do i=1, ubound(noprop,1)
        expcols=expcols+noprop(i)
      end do
      
      
            
      open(newunit=expfile, file=adjustl(trim(fileinputs)), iostat=ierr)
      
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
        do j=2, ubound(tmpdata,1)
          


      
      print *, expcols ; stop
      


        call file_error(expfile)
      end if
      
    
    end subroutine reader
    
    subroutine get_objval()
      use typy
      
      call reader()
    
    end subroutine get_objval
  


end module objfnc
