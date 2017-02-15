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

      class(pde_str), intent(in out) :: pde_loc
      integer :: i_err
      integer(kind=ikind) :: i, n
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      character(len=4096) :: msg

      
      

      pde_loc%problem_name(1) = "heat"
      pde_loc%problem_name(2) = "Heat conduction equation with convection (Sophocleous, 1979)"

      pde_loc%solution_name(1) = "temperature" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "T " !popisek grafu

      pde_loc%flux_name(1) = "heat_flux"  
      pde_loc%flux_name(2) = "heat flux [W.L-2]"

      pde_loc%mass_name(1) = "skip_this_file"
      pde_loc%mass_name(2) = "n/a"
      
      

      call find_unit(file_heat, 200)
      open(unit = file_heat, file="drutes.conf/heat/heat.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
	print *, "missing drutes.conf/heat/heat.conf file"
	ERROR STOP
      end if
     

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
			      

 
      write(unit=msg, fmt=*) "HINT 1: Are all values anisotropy defining anisotropical diffusivity positive? ", new_line("a"), &
	"HINT 2 : Have you defined enough values for anisotropy &
	(e.g. for 2D define angle and the maximal and minimal value of diffusivity, in total 3 values)?", new_line("a"),&
	"HINT 3: The number of lines with heat conductivity has to correspond to the number of materials & 
	defined by your mesh"
      
      
      allocate(tmp_array(drutes_config%dimen + 1))
      do i=1, ubound(heatpar,1)
	allocate(heatpar(i)%lambda_loc(drutes_config%dimen))
	call fileread(r=tmp_array, fileid=file_heat, ranges=(/0.0_rkind, huge(tmp)/), errmsg=trim(msg))
	heatpar(i)%anisoangle = tmp_array(1)
	heatpar(i)%lambda_loc = tmp_array(2:drutes_config%dimen + 1)
	allocate(heatpar(i)%lambda(drutes_config%dimen, drutes_config%dimen))
	call set_tensor(heatpar(i)%lambda_loc, (/heatpar(i)%anisoangle/), heatpar(i)%lambda)
      end do
      
      
      write(unit=msg, fmt=*) "Did you specify convection vector component for each coordinate (e.g. x,y,z)"
      do i=1, ubound(heatpar,1)
        allocate(heatpar(i)%convection(drutes_config%dimen))
        call fileread(r=heatpar(i)%convection, fileid=file_heat, errmsg=trim(msg))
      end do
        
      write(unit=msg, fmt=*) "Hint: The number of lines for the initial temperature has to be equal to the number of materials."
      do i=1, ubound(heatpar,1)
       call fileread(r=heatpar(i)%Tinit, fileid=file_heat, errmsg=trim(msg))
      end do
       
		    
      
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