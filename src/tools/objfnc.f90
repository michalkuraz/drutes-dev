module objfnc
  use typy
  
  private :: reader
  public :: get_objval
  
  integer, private :: exp_file
  
  real(kind=rkind), dimension(:,:), allocatable, private :: exp_data
  real(kind=rkind), dimension(:,:), allocatable, private :: model_data
  integer(kind=ikind), dimension(:), allocatable, private :: obs_ids, columns
  integer(kind=ikind), private :: noprop
  character(len=4096), private :: fileinputs
  
  contains
    
    subroutine reader()
      use typy
      use readtools
      use core_tools
      
      integer :: fileid
      
      integer(kind=ikind) :: n, i
      character(len=4096) :: msg
      
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
        call fileread(columns(i), fileid, ranges=(/0_ikind, huge(n)/), errmsg=msg)
      end do
      
      call fileread(fileinputs, fileid)
      
      call write_log("The file with inputs is", text2=trim(fileinputs))
      
    
    end subroutine reader
    
    subroutine get_objval()
      use typy
      
      call reader()
    
    end subroutine get_objval
  


end module objfnc
