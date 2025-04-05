
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

!> \file heat_reader.f90
!! \brief Reads configs for heat equation.
!<




module heat_reader
  public :: heat_read
  
  contains

    subroutine heat_read(pde_loc)
      use typy
      use globals
      use global_objs
      use core_tools
      use heat_globals
      use readtools
      use pde_objs
      use debug_tools

      class(pde_str), intent(in out) :: pde_loc
      integer :: i_err
      integer(kind=ikind) :: i, n, counter
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      character(len=4096) :: msg
      character(len=4) :: icondfile

      
      

      pde_loc%problem_name(1) = "heat"
      pde_loc%problem_name(2) = "Heat conduction equation with convection (Sophocleous, 1979)"

      pde_loc%solution_name(1) = "temperature" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "T " !popisek grafu

      pde_loc%flux_name(1) = "heat_flux"  
      pde_loc%flux_name(2) = "heat flux [W.L-2]"
      
      allocate(pde_loc%mass_name(0,2))

      call find_unit(file_heat, 200)
      open(unit = file_heat, file="drutes.conf/heat/heat.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
        print *, "missing drutes.conf/heat/heat.conf file"
        ERROR STOP
      end if
     
      call fileread(with_richards, file_heat)

      allocate(heatpar(maxval(elements%material)))
      
      call fileread(n, file_heat)
      
      backspace(file_heat)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/heat/heat.conf  &
        the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_heat, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
        errmsg=trim(msg))
	
      write(unit=msg, fmt=*) "HINT 1: Is the heat capacity (matrix/matrix) positive?", new_line("a"), &
        "   HINT 2 : Is the number of heat capacity values corresponding to the amount of layers?"
	
      do i=1, ubound(heatpar,1)
        call fileread(heatpar(i)%C, file_heat, ranges=(/0.0_rkind,huge(0.0_rkind)/),&
        errmsg=trim(msg))
      end do
      
      write(unit=msg, fmt=*) "HINT 1: Is the heat capacity (matrix/water) positive?", new_line("a"), &
        "   HINT 2 : Is the number of heat capacity values corresponding to the amount of layers?"
        
      do i=1, ubound(heatpar,1)
        call fileread(heatpar(i)%C_w, file_heat, ranges=(/0.0_rkind,huge(0.0_rkind)/),&
        errmsg=trim(msg))
      end do
      
      
      write(unit=msg, fmt=*) "HINT 1: Is the heat conductivity positive?", new_line("a"), &
        "   HINT 2 : Is the number of heat conductivity values corresponding to the amount of layers?"
                  

 

      write(unit=msg, fmt=*) "HINT 1: Are all values anisotropy defining anisotropical diffusivity positive? ", &
       new_line("a"), new_line("a"),  &
        "   HINT 2: Have you defined an EXACT NUMBER of values for anisotropy?", new_line("a"), &
        "     (e.g. for 2D define angle and the maximal and minimal value of diffusivity, in total 3 values", new_line("a"),&
         "      if you switch from 2D to 1D don't forget to erase the last value, otherwise this error is generated)", &
        new_line("a"), new_line("a"), &
        "   HINT 3: The number of lines with heat conductivity has to correspond to the number of materials & 
        defined by your mesh", new_line("a"), new_line("a")

      
      
      allocate(tmp_array(drutes_config%dimen + 1))
      do i=1, ubound(heatpar,1)
        allocate(heatpar(i)%lambda_loc(drutes_config%dimen))

        call fileread(r=tmp_array, fileid=file_heat, ranges=(/0.0_rkind, huge(tmp)/), errmsg=trim(msg), checklen=.TRUE.)

        heatpar(i)%anisoangle = tmp_array(1)
        heatpar(i)%lambda_loc = tmp_array(2:drutes_config%dimen + 1)
        allocate(heatpar(i)%lambda(drutes_config%dimen, drutes_config%dimen))
        call set_tensor(heatpar(i)%lambda_loc, (/heatpar(i)%anisoangle/), heatpar(i)%lambda)
      end do
      
      
      if (.not. with_richards) then
        write(unit=msg, fmt=*) "Did you specify convection vector component for each coordinate (e.g. x,y,z)"
        do i=1, ubound(heatpar,1)
          allocate(heatpar(i)%convection(drutes_config%dimen))
          call fileread(r=heatpar(i)%convection, fileid=file_heat, errmsg=trim(msg))
        end do
      end if
        
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      do i=1, ubound(heatpar,1)
        call fileread(icondfile, fileid=file_heat)
        if (icondfile /= "file") then
          backspace file_heat
          call fileread(r=heatpar(i)%Tinit, fileid=file_heat, errmsg=trim(msg))
        else
          Tinitfile = .true.
          EXIT
        end if
      end do
      if (Tinitfile) then
        open(newunit=file_heaticond, file="drutes.conf/heat/heaticond1D.in", action="read", status="old", iostat=i_err)
        if (i_err /= 0) then
          write(unit=msg, fmt=*) "Unable to open file drutes.conf/heat/heaticond1D.in", new_line("a"), &
                                 "if you don't want to use this file for heat initial condition definition", new_line("a"), &
                                 "change value in drutes.conf/heat/heat.conf", new_line("a"), &
                                 "instead of a keyword file for initial condition set some concrete value",  new_line("a"), &
                                 "DRUtES will then use piece-wise constant distribution of the initial condition."
          call file_error(file_heat, errmsg = trim(msg))
        end if
        counter = 0
        do 
          call comment(file_heaticond)
          read(unit=file_heaticond, fmt=*, iostat=i_err) tmp, tmp
          if (i_err == 0) then
            counter = counter + 1
          else
            EXIT
          end if
        end do
        if (counter < 2) then
          write(unit=msg, fmt=*) "You have specified less than two points for interpolation."
          call file_error(file_heaticond, errmsg=trim(msg))
        end if
        
        close(file_heaticond)
        
        open(newunit=file_heaticond, file="drutes.conf/heat/heaticond1D.in", action="read", status="old")
        
        allocate(Tinitpoints(counter,2))
        
        write(unit=msg, fmt=*) "Have you specified x and T values (two columns)?"
        
        do i=1, counter
          call fileread(Tinitpoints(i,:), file_heaticond, errmsg=trim(msg))
        end do
          
      end if
       
		    
      
      write(unit=msg, fmt=*) "Hint: The number of lines for the heat source has to be equal to the number of materials."
      do i=1, ubound(heatpar,1)
       call fileread(r=heatpar(i)%source, fileid=file_heat, errmsg=trim(msg))
      end do
      

      write(unit=msg, fmt=*) "The number of boundaries should be greater than zero and smaller or equal the number of nodes"
      
      call fileread(n, file_heat, ranges=(/1_ikind, nodes%kolik/),&
        errmsg=trim(msg))

      call readbcvals(unitW=file_heat, struct=pde_loc%bc, dimen=n, &
          dirname="drutes.conf/heat/")
          
          

      
      
      

    end subroutine heat_read		
    



  

end module heat_reader
