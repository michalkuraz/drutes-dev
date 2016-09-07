! Copyright 2008 Michal Kuraz, Petr Mayer

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

!> \mainpage 
!! \section  Introduction (DRUtES = Dual Richards' Unsaturated Equation Solver)
!!  Unsaturated flow (moisture movement) equation solver with preferential flow support and convection-dispersion-reaction equation.
!<



program main
  use typy
  use drutes_init
  use globals
  use core_tools
  use manage_pointers
  use fem
  use postpro
  use decomposer
  use read_inputs
  use feminittools
  use debug_tools
  use simplelinalg
  use pde_objs
  use debug_tools

  character(len=256) :: writer
  character(len=2)   :: ch
  logical :: success
  real ::  stop_time
  real(kind=rkind) :: r
  integer :: fileid, i

  
  
  call system("rm -rf out/*")
  
  
  if (this_image() == 1) then
  
    call getcwd(dir_name)

    call cpu_time(start_time)



    version_id%number = "1.201601/"
    version_id%reliability = "beta "
    
    call get_cmd_options()
    
    terminal = 6
    
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)"This program is free software: you can redistribute it and/or modify"
    write(unit=terminal, fmt=*)"it under the terms of the GNU General Public License as published by"
    write(unit=terminal, fmt=*)"the Free Software Foundation, either version 3 of the License, or"
    write(unit=terminal, fmt=*)"(at your option) any later version."
    write(unit=terminal, fmt=*)"This program is distributed in the hope that it will be useful,"
    write(unit=terminal, fmt=*)"but WITHOUT ANY WARRANTY; without even the implied warranty of"
    write(unit=terminal, fmt=*)"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the"
    write(unit=terminal, fmt=*)"GNU General Public License for more details."
    write(unit=terminal, fmt=*)"You should have received a copy of the GNU General Public License"
    write(unit=terminal, fmt=*)"along with this program. If not, see <http://www.gnu.org/licenses/>."
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)" "
    write(unit=terminal, fmt=*)" "
    
    write(unit=terminal, fmt=*) " " //achar(27)//'[94m', "DRUtES" //achar(27)//'[0m', &
	   " version: " //achar(27)//'[92m', version_id, " " //achar(27)//'[0m'
	   
    write(unit=terminal, fmt=*)" "
    write(unit=terminal, fmt=*)" " 
    
    print *, " " //achar(27)//'[94m', "DRUtES" //achar(27)//'[0m', &
             " version: " //achar(27)//'[92m', version_id, " " //achar(27)//'[0m'
      
    call parse_globals() 
    
    call init_measured()
        
    call write_log("number of nodes:", int1=nodes%kolik, text2="number of elements:", int2=elements%kolik)
        
    call set_pointers()
    
    call init_observe()
        
    call feminit()
    
        
    if (drutes_config%it_method == 1 .or. drutes_config%it_method == 2) then
      call init_decomp()
    end if

  end if
  
  call write_log("DRUtES solves ", text2=adjustl(trim(drutes_config%fullname)))

  call solve_pde(success)      
  

  sync all
  
  if (this_image() == 1) then
    call cpu_time(stop_time)
    
    
    select case (int(stop_time - start_time))
            case (0:60)	
                    call write_log(text="# real elapsed CPU time =", real1=1.0_rkind*(stop_time - start_time), text2="s")
            case (61:3600)
		    call write_log(text="# real elapsed CPU time =", real1=(stop_time - start_time)/60.0_rkind, text2="min")
            case (3601:86400)
		    call write_log(text="# real elapsed CPU time =", real1=(stop_time - start_time)/3600.0_rkind, text2="hrs")
            case default
		    call write_log(text="# real elapsed CPU time =", real1=(stop_time - start_time)/86400.0_rkind, text2="days") 
   end select
   
    call find_unit(fileid)
    open(unit=fileid, file="out/cpu.time", action="write", status="replace")
    write(unit=fileid, fmt=*) (stop_time - start_time)
    close(fileid)

    call write_log("minimal adjusted time step was", real1=minimal_dt)

    call write_log("F I N I S H E D !!!!!")
    
    
    print *, "DRUTES version: ", version_id, "terminated normally"
  end if

   sync all

   do i=1, NUM_IMAGES()
    call get_RAM_use()
    call flush(6)
   end do
   

end program main
