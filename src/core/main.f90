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

!> \file main.f90
!! \brief Main program
!<

!> \mainpage 
!! \section DRUtES
!! originally developed as Dual Richards' Unsaturated Equation Solver
!! currently DRUtES is an object oriented Fortran 2003/2008 library for solving coupled nonlinear convection-diffusion reaction equation in the following format: \n
!! \f[ \begin{split} \sum_{i=1}^n C_{1,i} \frac{\partial p_i}{\partial t} &= \sum_{i=1}^n \left( \nabla \cdot \mathbf{D}_{1,i} \nabla p_i - \nabla  \cdot (\vec{q}_{1,i} p_i) -  \sum_{r=0}^{r_{max}} \lambda_{1,i}p_i^r \right)  \\ \\ & \vdots  \\  \sum_{i=1}^n C_{n,i} \frac{\partial p_i}{\partial t} &= \sum_{i=1}^n \left( \nabla \cdot \mathbf{D}_{n,i} \nabla p_i - \nabla  \cdot (\vec{q}_{n,i} p_i) -  \sum_{r=0}^{r_{max}} \lambda_{n,i}p_i^r \right) \end{split}   \f]
!! \n
!! where \f$n\f$ is allocatable dimension of the PDE problem,  \f$r_{max}\f$ is an arbitrary reaction order, and  \f$p_i\f$ is the  \f$i\f$ component of the solution of the PDE problem.
!!  
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
  use re_analytical
  use objfnc
  use printtools
  
  character(len=256) :: writer
  character(len=2)   :: ch
  logical :: success
  real ::  stop_time
  real(kind=rkind) :: r, t
  integer :: fileid, i, j
  
  
  call system("rm -rf out/*")
  
  
  if (this_image() == 1) then
  
    call getcwd(dir_name)

    call cpu_time(start_time)



    version_id%number = "8.0/2018"
    version_id%reliability = "beta "
    
    call get_cmd_options()
    
    terminal = 6
    
    call print_logo()
    
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)"DRUtES is a free software: you can redistribute it and/or modify"
    write(unit=terminal, fmt=*)"it under the terms of the GNU General Public License as published by"
    write(unit=terminal, fmt=*)"the Free Software Foundation, either version 3 of the License, or"
    write(unit=terminal, fmt=*)"(at your option) any later version."
    write(unit=terminal, fmt=*)"DRUtES is distributed in the hope that it will be useful,"
    write(unit=terminal, fmt=*)"but WITHOUT ANY WARRANTY; without even the implied warranty of"
    write(unit=terminal, fmt=*)"MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the"
    write(unit=terminal, fmt=*)"GNU General Public License for more details."
    write(unit=terminal, fmt=*)"See <http://www.gnu.org/licenses/>."
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)"---------------------------------------------------------------------------"
    write(unit=terminal, fmt=*)" "
    write(unit=terminal, fmt=*)" "
    
    write(unit=terminal, fmt=*) " " //achar(27)//'[94m', "DRUtES" //achar(27)//'[0m', &
	   " version: " //achar(27)//'[92m', version_id, " " //achar(27)//'[0m'
	   
    write(unit=terminal, fmt=*)" "
    write(unit=terminal, fmt=*)" " 
    
      
       
    call parse_globals() 
    
    call init_measured()
  
    call write_log("number of nodes:", int1=nodes%kolik, text2="number of elements:", int2=elements%kolik)
        
    call set_pointers()
 
    call init_observe()
   
    call feminit()
    
    if (drutes_config%it_method == 1 .or. drutes_config%it_method == 2) then
      call init_decomp()
    end if
    
     if (objval%compute) call objval%read_config()

  end if
  
  call write_log("DRUtES solves ", text2=adjustl(trim(drutes_config%fullname)))

  call solve_pde(success)    

  
  if (objval%compute) call objval%getval()
  
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
