module printtools
  public :: progressbar
  public :: time2finish
  public :: printtime

  contains
  !> if default output is stdout, then this procedure prints nice progress bar

    subroutine progressbar(j)
      use globals
  
      integer, intent(in) :: j
!       character(len=*), intent(in), optional :: text
      integer, save :: print_pos
      integer :: k  
      character(len=62):: bar

      if (j < 1) then
	print_pos = -1
      end if


      if (terminal == 6) then
	if (j > print_pos) then
	  bar="???% |                                                  |"  
	  write(unit=bar(1:3),fmt="(i3)") j  
	  do k=1, nint(j/2.0)  
	    bar(6+k:6+k)="*"  
	  end do  
	  ! print the progress bar.  
! 	  if (present(text) ) then
! 	    write(unit=terminal,fmt="(a, a,a1,a63)", advance="no") "  ", "neco", "+",char(13), bar 
! 	  else
	    write(unit=terminal,fmt="(a1,a1,a63)", advance="no") "+",char(13), bar  
! 	  end if
	    
	  call flush(terminal)
	  print_pos = print_pos + 1
	  if (j>= 100) then
	    write(unit=terminal, fmt=*) "  "
	    print_pos = 0
	  end if
	end if
	return
      else
	return
      end if

    end subroutine progressbar
    
    subroutine printtime(message, t)
      use globals
      character(len=*), intent(in) :: message
      real(4), intent(in) :: t
      
      select case (int(t))
	case (0:60)   
	  write(unit=terminal, fmt="(a, F10.2, a)") trim(message), t, " s"
	case (61:3600)
	  write(unit=terminal, fmt="(a, F10.2, a)") trim(message), t/60, " min"
	case (3601:86400)
	  write(unit=terminal, fmt="(a, F10.2, a)") trim(message), t/3600, " hrs"
	case (86401:2592000)
	  write(unit=terminal, fmt="(a, F10.2, a)") trim(message), t/86400, " days"
	case (2592001:31536000)	
	  write(unit=terminal, fmt="(a, F10.2, a)") trim(message), t/2592000, " months"
	case default
	  write(unit=terminal, fmt=* ) trim(message), t/31536000, " " //achar(27)//'[91m',  " YEARS !!!!...:)" &
						//achar(27)//'[0m'		    
      end select
      
    end subroutine printtime


    subroutine time2finish()
      use typy
      use globals
      use core_tools

      real(kind=rkind) :: tmp, tmp1
      real :: konec
      integer :: fileid

      call cpu_time(konec)
      
      tmp1 = min(100.0,time/end_time*100)
      

      tmp = max(0.0,(end_time-time)/(time/(konec-start_time)))

      select case (int(tmp))
	case (0:60)   
	  write(unit=terminal, fmt="(a, F10.2, a)") " #time to finish =", tmp, " s"
	case (61:3600)
	  write(unit=terminal, fmt="(a, F10.2, a)") " #time to finish =", tmp/60, " min"
	case (3601:86400)
	  write(unit=terminal, fmt="(a, F10.2, a)") " #time to finish =", tmp/3600, " hrs"
	case (86401:2592000)
	  write(unit=terminal, fmt="(a, F10.2, a)") " #time to finish =", tmp/86400, " days"
	case (2592001:31536000)	
	  write(unit=terminal, fmt="(a, F10.2, a)") " #time to finish =", tmp/2592000, " months"
	case default
	  write(unit=terminal, fmt=* ) " #time to finish =", tmp/31536000, " " //achar(27)//'[91m',  " YEARS !!!!...:)" &
						//achar(27)//'[0m'		    
      end select
      write(unit=terminal, fmt="(F10.2, a)") tmp1, "% done"
      
      if (www) then
	call find_unit(fileid,1000)
	open(unit=fileid, file="4www/progress", status="replace", action="write")
	write(unit=fileid, fmt=*) tmp1
	select case (int(tmp))
	  case (0:60)   
	    write(unit=fileid, fmt="(a, F10.2, a)") " #time to finish =", tmp, " s"
	  case (61:3600)
	    write(unit=fileid, fmt="(a, F10.2, a)") " #time to finish =", tmp/60, " min"
	  case (3601:86400)
	    write(unit=fileid, fmt="(a, F10.2, a)") " #time to finish =", tmp/3600, " hrs"
	  case (86401:2592000)
	    write(unit=fileid, fmt="(a, F10.2, a)") " #time to finish =", tmp/86400, " days"
	  case (2592001:31536000)	
	    write(unit=fileid, fmt="(a, F10.2, a)") " #time to finish =", tmp/2592000, " months"
	  case default
	    write(unit=fileid, fmt=* ) " #time to finish =", tmp/31536000, " " //achar(27)//'[91m',  " YEARS !!!!...:)" &
						  //achar(27)//'[0m'		    
	end select
	close(fileid)
      end if
 

    end subroutine time2finish

end module printtools