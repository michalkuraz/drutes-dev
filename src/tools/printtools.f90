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

!> \file printtools.f90
!! \brief Terminal print functions.
!<



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
      use debug_tools

      real(kind=rkind) :: tmp, tmp1
      real :: now
      integer :: fileid, ierrtime
      
      
      open(newunit=fileid, file="/proc/uptime", action="read", status="old", iostat=ierrtime)
    
    
      if (ierrtime /= 0) then
        call cpu_time(now)
      else
        read(unit=fileid, fmt=*) now
        close(fileid)
      end if

      
      tmp1 = min(100.0,time/end_time*100)
      
      tmp = max(0.0,(end_time-time)/(time/(now-start_time)))

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
        select case (nint(tmp))
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
    
    subroutine print_logo(fileid)
      use globals
      
      integer, intent(in), optional :: fileid
      
      integer :: iunit
      
      if (present(fileid)) then
        iunit = fileid
      else
        iunit=terminal
      end if
      
      write(unit=iunit, fmt=*) "#              _____________________  _______________________             #"
      write(unit=iunit, fmt=*) "#              ___  __ \__  __ \_  / / /_  /___  ____/_  ___/             #"
      write(unit=iunit, fmt=*) "#              __  / / /_  /_/ /  / / /_  __/_  __/  _____ \              #"
      write(unit=iunit, fmt=*) "#              _  /_/ /_  _, _// /_/ / / /_ _  /___  ____/ /              #"
      write(unit=iunit, fmt=*) "#              /_____/ /_/ |_| \____/  \__/ /_____/  /____/               #"
      write(unit=iunit, fmt=*) "#                                                                         #"
      write(unit=iunit, fmt=*) " "
      
  end subroutine print_logo
      

end module printtools
