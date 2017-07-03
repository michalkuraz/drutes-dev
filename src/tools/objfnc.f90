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
  integer(kind=ikind), dimension(:), allocatable, private :: obs_ids, columns
  integer(kind=ikind), private :: noprop
  character(len=4096), private :: fileinputs
  type(ram_limit_str), private :: ram_limit
  
  
  
  contains
    
    subroutine reader()
      use typy
      use readtools
      use core_tools
      
      integer :: fileid, expfile, ierr
      
      integer(kind=ikind) :: n, i, counter
      character(len=4096) :: msg
      real(kind=rkind) :: r
      real(kind=rkind), dimension(:), allocatable :: tmpdata
      
      open(newunit=fileid, file="drutes.conf/inverse_modeling/objfnc.conf")
      
      call fileread(n, fileid, ranges=(/0_ikind, huge(n)/))
      
      allocate(obs_ids(n))
      
      msg="Is the number of observation points for evaluating your objective function correct?"
      
      do i=1, n
        call fileread(obs_ids(i), fileid, ranges=(/0_ikind, huge(n)/), errmsg=msg)
      end do
        
      call fileread(noprop, fileid, ranges=(/0_ikind, huge(n)/))
      
      allocate(columns(noprop))
      
      msg="Is the number of properties for evaluating your objective function correct?"
      
      do i=1, noprop
!       time, val, massval, advectval(1:D), observation_array(i)%cumflux(proc) - in total 4 properties + time
        call fileread(columns(i), fileid, ranges=(/0_ikind, 4_ikind/), errmsg=msg)
      end do
      
      call fileread(fileinputs, fileid)
      
      call write_log("The file with inputs for inverse modeling is", text2=trim(fileinputs))

      call fileread(ram_limit%set, fileid)
      
      call fileread(ram_limit%memsize, fileid)
      
      call fileread(ram_limit%units, fileid, options=(/"kB", "MB", "GB"/))
      
      open(newunit=expfile, file=adjustl(trim(fileinputs)), iostat=ierr)
      
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
        read(unit=expfile, iostat=ierr) r
        if (ierr /= 0) then
          counter = counter-1
          EXIT
        end if
      end do
      
      allocate(exp_data(noprop))
      
      allocate(exp_data(1)%time(counter))
      
      do i=1, noprop
        allocate(exp_data(i)%data(counter))
      end do
      
      close(exp_file)
      
      open(newunit=expfile, file=adjustl(trim(fileinputs)))

!       do i=1, counter
!         call fileread( tmpdata(, exp
!       
      
      
    
    end subroutine reader
    
    subroutine get_objval()
      use typy
      
      call reader()
    
    end subroutine get_objval
  


end module objfnc
