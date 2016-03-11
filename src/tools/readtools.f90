module readtools
  public :: file_error
  public :: fileread
  private :: read_int, read_int_std, read_int_array, read_real,  read_real_array, read_char
  public :: comment
  public :: readbcvals
  public :: set_tensor

  interface fileread
    module procedure read_int
    module procedure read_int_std
    module procedure read_int_array
    module procedure read_real
    module procedure read_real_array
    module procedure read_char
    module procedure read_logical
  end interface fileread

  contains

    subroutine read_int(i, fileid, errmsg, ranges, minimal, maximal, noexit)
      use typy
      use globals
      integer(kind=ikind), intent(out) :: i
      integer, intent(in) :: fileid
      character(len=*), intent(in), optional :: errmsg
      integer(kind=ikind), dimension(:), intent(in),  optional :: ranges
      integer, intent(in), optional :: minimal, maximal
      logical, intent(in), optional :: noexit

      !local vars
      integer :: ierr

      call comment(fileid)

      read(unit=fileid, fmt=*, iostat=ierr) i

      if (ierr /= 0) then
        if (present(errmsg)) then
          call file_error(fileid, errmsg)
        end if
        call file_error(fileid)
      end if

      if (present(ranges)) then
        if (i<ranges(1) .or. i> ranges(2)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
          if (.not. present(noexit) .or. .not. noexit) then
	    if (present(errmsg)) then
	      call file_error(fileid,errmsg)
	    else
	      call file_error(fileid)
	    end if
	  end if
        end if
      end if

    end subroutine read_int


    subroutine read_int_std(i, fileid, errmsg, ranges, noexit)
      use typy
      use globals
      integer(kind=lowikind), intent(out) :: i
      integer, intent(in) :: fileid
      character(len=*), intent(in), optional :: errmsg
      integer(kind=ikind), dimension(:), intent(in),  optional :: ranges
      logical, intent(in), optional :: noexit


      !local vars
      integer :: ierr

      call comment(fileid)

      read(unit=fileid, fmt=*, iostat=ierr) i

      if (ierr /= 0) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
        if (.not. present(noexit) .or. .not. noexit) then
	  if (present(errmsg)) then
	    call file_error(fileid,errmsg)
	  else
	    call file_error(fileid)
	  end if
	end if
      end if
      
      
       if (present(ranges)) then
        if (i<ranges(1) .or. i> ranges(2)) then
          write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
          if (.not. present(noexit) .or. .not. noexit) then
	    if (present(errmsg)) then
	      call file_error(fileid,errmsg)
	    else
	      call file_error(fileid)
	    end if
	  end if
        end if
      end if

    end subroutine read_int_std


    subroutine read_int_array(i, fileid, errmsg, noexit)
      use typy
      integer(kind=ikind), dimension(:), intent(out) :: i
      integer, intent(in) :: fileid
      character(len=*), intent(in), optional :: errmsg
      logical, intent(in), optional :: noexit

      !local vars
      integer :: ierr

      call comment(fileid)

      read(unit=fileid, fmt=*, iostat=ierr) i

      if (ierr /= 0) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
        if (.not. present(noexit) .or. .not. noexit) then
	  if (present(errmsg)) then
	    call file_error(fileid,errmsg)
	  else
	    call file_error(fileid)
	  end if
	end if
      end if

    end subroutine read_int_array
    
    
    subroutine read_real(r, fileid, ranges, errmsg, noexit)
      use typy
      use globals
      real(kind=rkind), intent(out) :: r
      integer, intent(in) :: fileid
      real(kind=rkind), dimension(:), optional :: ranges
      character(len=*), intent(in), optional :: errmsg
      logical, intent(in), optional :: noexit
      
      !logical vars
      integer :: ierr
      
      call comment(fileid)
      
      read(unit=fileid, fmt=*, iostat=ierr) r
      
      if (present(ranges)) then
	if (r < ranges(1) .or. r > ranges(2)) then
	  write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
	  ierr = -1
	end if
      end if
      
      if (ierr /= 0) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
        if (.not. present(noexit) .or. .not. noexit) then
	  if (present(errmsg)) then
	    call file_error(fileid,errmsg)
	  else
	    call file_error(fileid)
	  end if
	end if
      end if

    end subroutine read_real
    
    

    subroutine read_real_array(r, fileid, ranges, errmsg, noexit)
      use typy
      use globals
      real(kind=rkind), dimension(:), intent(out) :: r
      integer, intent(in) :: fileid
      real(kind=rkind), dimension(:), optional :: ranges
      character(len=*), intent(in), optional :: errmsg
      logical, intent(in), optional :: noexit
      
            !logical vars
      integer :: ierr
      integer(kind=ikind) :: i
      
      call comment(fileid)
      
      read(unit=fileid, fmt=*, iostat=ierr) r
      
      if (ierr /= 0) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
        if (.not. present(noexit) .or. .not. noexit) then
	  if (present(errmsg)) then
	    call file_error(fileid,errmsg)
	  else
	    call file_error(fileid)
	  end if
	end if
      end if
      
      if (present(ranges)) then
	do i=1, ubound(r,1)
	  if (r(i) < ranges(1) .or. r(i) > ranges(2)) then
	    write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "incorrect ranges of input parameter(s)" //achar(27)//'[0m'
	    call file_error(fileid, errmsg)
	  end if
	end do
      end if
      
      
    end subroutine read_real_array
    
    
    
    subroutine read_char(ch, fileid, errmsg, options, noexit)
      character(len=*), intent(out) :: ch
      integer, intent(in) :: fileid
      character(len=*), intent(in), optional :: errmsg
      character(len=*), dimension(:), intent(in), optional :: options
      logical, intent(in), optional :: noexit

      integer :: ierr, i
      logical :: ok
      
      call comment(fileid)
      
      read(unit=fileid, fmt=*, iostat=ierr) ch
      
      if (present(options)) then
	ok = .false.
	do i=1, ubound(options, 1)
	  if (adjustl(trim(ch)) == adjustl(trim(options(i)))) then
	    ok = .true.
	    EXIT
	  end if
	end do
      else
	ok = .true.	
      end if
      
      if (ierr /= 0 .or. .not.(ok) .and. (.not. present(noexit) .or. .not. noexit)) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
	if (present(errmsg)) then
	  call file_error(fileid,errmsg)
	else
	  call file_error(fileid)
	end if
      end if
      
      
      
      
    end subroutine read_char
    
    subroutine read_logical(l, fileid, errmsg, noexit)
      use globals
      logical, intent(out) :: l
      integer, intent(in) :: fileid
      character(len=*), intent(in), optional :: errmsg
      logical, intent(in), optional :: noexit

      integer :: ierr
      character(len=1) :: ch
      
      call comment(fileid)
      
      read(unit=fileid, fmt=*, iostat=ierr) ch

      if (ierr == 0) then 
	select case(ch)
	  case("y")
	    l = .true.
	  case("n")
	    l = .false.
	  case default
	    ierr = -1
	    write(unit=terminal, fmt=*) " " //achar(27)//'[91m', "You have defined incorrect value for logical parameter)" &
	     //achar(27)//'[0m'
	 end select
       end if
      
      if (ierr /= 0 .and. (.not. present(noexit) .or. .not. noexit)) then
        if (present(errmsg)) then
          call file_error(fileid,errmsg)
        end if
	if (present(errmsg)) then
	  call file_error(fileid,errmsg)
	else
	  call file_error(fileid)
	end if
      end if
      
      
    end subroutine read_logical
      
      

    subroutine file_error(iunit, message)
      use core_tools
      use globals

      integer, intent(in) :: iunit
      character(len=*), intent(in), optional :: message

      character(len=512) :: filename
      character(len=256) :: iaction, numero
      integer :: line, i, i_err, i3, i4, i5, fileid, err_read
      integer(4) :: i1, i2
      character(len=4096) :: string 

      string = "last value"
      i1=-100
      i2=-100
      inquire(unit=iunit, name=filename, action=iaction) 
      line = 1
      backspace(unit=iunit, iostat = i_err)

      read(unit=iunit, fmt=*, iostat = err_read) string
      
      
      do
        i1 = ftell(iunit)
        backspace (unit=iunit, iostat=i_err)
        line = line + 1
        i2 = ftell(iunit)
        if (i1 == i2) then
          EXIT
        end if
      end do

      if (line <= 2) then
        print *, "-------------------------------------------------------------------"
        print *, "WARNING: the line number was identified incorrectly, "
        print *, "likely it is a compiler error in ftell function, "
        print *, "find the input error on your own. :("
      end if

      print *, "-------------------------------------------------------------------"

      print *, "  "
      print *, " " //achar(27)//'[91m', "!!ERROR!!" //achar(27)//'[0m', &
                  "  in reading from file: ", trim(filename), " near line: ", line-1, &
                      " opened with status: ", trim(iaction),  " -> exiting DRUtES"
      print *, "  "
      

!       if (i_err == 0) then 
	if (err_read == 0) write(unit=terminal, fmt=*) "the value you have typed is: ", trim(string)
	write(unit=terminal, fmt=*) "--------------------------------"
!       end if
      
      
      if (present(message)) then
	write(unit=terminal, fmt=*) " " //achar(27)//'[91m', trim(message) //achar(27)//'[0m'
      end if
	
      
      ERROR STOP


    end subroutine file_error

    subroutine comment(Unit)
          
      integer, intent(in) :: Unit
      !character(len=1), intent(in), optional :: mark
      character(len=1) :: symbol
      character(len=1) :: String
      integer :: i_err
      
          
      !if (present(mark)) then
      !   symbol = mark
      !else
        symbol = "#"
    ! end if
        
      do
        read(unit=Unit,fmt = *, iostat = i_err ) String
        if (i_err < 0) then
          RETURN
        end if
        if (String == symbol) then
          CONTINUE
        else
          backspace Unit
          RETURN
        endif
      end do
          
    end subroutine comment


    !> this procedure set up an anisothropy tensor, the local anisothrophy values are supllied in values array, angle is the angle between the local and global system, tensor is the resulting secodn order tensor
    subroutine set_tensor(values, angle, tensor)
      use typy
      use globals
      use global_objs
      !>local anisothoropy values with respect to local x, local y and local z
      real(kind=rkind), dimension(:), intent(in) :: values
      !> angle between the local and global system of coordinates
      real(kind=rkind), dimension(:), intent(in) :: angle
      !> resulting second order tensor with respect to the global system of coordinates
      real(kind=rkind), dimension(:,:), intent(out) :: tensor
      real(kind=rkind), dimension(:,:), allocatable :: T
      integer(kind=ikind) :: i,j


      allocate(T(drutes_config%dimen, drutes_config%dimen))


      select case(drutes_config%dimen)
        case(1)
          tensor(1,1) = values(1)
          RETURN
        case(2)

          T(1,1) = cos(4.0_dprec*atan(1.0_dprec)/180.0_rkind*angle(1))
          T(1,2) = cos(4.0_dprec*atan(1.0_dprec)/180.0_rkind*(90-angle(1)))
          T(2,1) = cos(4.0_dprec*atan(1.0_dprec)/180.0_rkind*(90+angle(1)))
          T(2,2) = cos(4.0_dprec*atan(1.0_dprec)/180.0_rkind*(angle(1)))

        case(3)
          print *, "tensors not implemented in 3D"
          ERROR STOP
      end select            

      tensor = 0
      do i=1, ubound(tensor,1)
        tensor(i,i) = values(i)
      end do

      tensor = matmul(matmul(T, tensor), transpose(T))

      deallocate(T)

    end subroutine  set_tensor


    subroutine readbcvals(unitW, struct, dimen, dirname)
      use typy
      use globals
      use global_objs
      use pde_objs
      use core_tools

      integer, intent(in) :: unitW
      type(boundary_vals), dimension(:), allocatable, intent(out) :: struct
      integer(kind=ikind), intent(in) :: dimen
      !> directory name, where data with boundary condition are stored
      character(len=*), intent(in) :: dirname
      integer(kind=ikind) :: i, j
      character(len=3) :: ch
      character(len=1) :: y_file, y_cond
      integer :: ierr, fileid
      integer(kind=ikind) :: counter
      character(len=256) :: filename
      real(kind=rkind) :: tmp
      character(len=1024) :: msg
      
      
      if (.not. allocated(struct)) then
	allocate(struct(101:101+dimen-1))
      else
	print *, "W: strange bc struct already allocated with bounds:", lbound(struct,1), ":", ubound(struct,1)
      end if
      
      if (maxval(nodes%edge) /= ubound(struct,1)) then
	inquire(unit=unitW, name=filename)
	if (maxval(nodes%edge)-100 - dimen == ubound(measured_pts,1) ) then
	  write(unit=msg, fmt="(a)") "There is an inconsistent boundary description, it seems like you have forgot & 
	       to include the points with measurement into boundary definition"
	else if (maxval(nodes%edge) < 100) then
	  write(unit=msg, fmt=*) "Your mesh does not define any boundaries, you are completely wrong!!"
	else 			!a,I3,a,I3,a,I3,a,I3,a,I3,a
	  write(unit=msg, fmt="(a,I3,a,I3,a,I3,a,a,a,I3,a)")  "There is an inconsistent boundary description, &
								  the mesh configuration defines "& 
	  , maxval(nodes%edge)-100 - ubound(measured_pts,1) , " boundary(ies), your measurement points define an extra", &
	  ubound(measured_pts,1), " boundary(ies), thus in total you need:", &
	  maxval(nodes%edge)-100, " boundary(ies), however the input model file: ", &
	  trim(filename), " defines ", dimen, " boundary(ies)."

	end if
	call file_error(unitW, trim(msg))
      end if

      
      do i = lbound(struct,1), ubound(struct,1)
        call comment(unitW)
        read(unit=unitW, fmt=*, iostat=ierr) struct(i)%ID, struct(i)%code, y_file, struct(i)%value, struct(i)%layer!, y_cond
        select case(y_file)
          case("y")
            struct(i)%file = .true.
          case("n")
            struct(i)%file = .false.
          case default
            call file_error(unitW, &
            "set [y/n] value for file in boundary condition description, see file info bellow")
        end select


            struct(i)%icond4neumann = .false.

        if (ierr /= 0) then
          inquire(unit=unitW, name=filename)
          write(unit=msg, fmt=*) "HINT: check number of boundary records in file: ", trim(filename)
          call file_error(unitW, trim(msg))
        end if

        if (struct(i)%file) then
          call find_unit(fileid)
          write(unit=filename, fmt="(a, I3, a)") trim(dirname), i, ".bc"
          open(unit=fileid, file=trim(filename), action="read", status="old", iostat=ierr)
          if (ierr /= 0) then
            write(unit=msg, fmt=*) "ERROR: if using unsteady boundary value data, you must supply the file with data!", &
            "ERROR: file required:", trim(filename)
            call file_error(unitW)
          end if
          counter = 0
          do 
            counter = counter + 1
            read(unit=fileid, fmt=*, iostat=ierr) tmp
            if (ierr /= 0) then
              EXIT
            end if
          end do

          do j=1, counter+1
            backspace fileid
          end do

          counter = counter - 1
          
          allocate(struct(i)%series(counter,2))

          do j=1, counter
            read(unit=fileid, fmt=*, iostat=ierr) struct(i)%series(j,:)
            if (ierr /= 0) then
              write(unit=msg, fmt="(a,a,a,i4)"), "incorrect data in file:  ", trim(filename), "    at line:", j
              call file_error(fileid, msg)
            end if
          end do

          struct(i)%series_pos = 1


          if (struct(i)%series(counter,1) < end_time) then
            call file_error(fileid, &
            "ERROR: end time in unsteady boundary data must be equal or higher than the simulation end time" )
          end if


          close(fileid)
        end if

      end do

      

    end subroutine readbcvals

end module readtools