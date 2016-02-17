module ADE_reader
  public :: ADE_read
  
  contains

    subroutine ADE_read()
      use typy
      use globals
      use global_objs
      use core_tools
      use ADE_globals
      use readtools
      use pde_objs

      integer :: i_err
      integer(kind=ikind) :: i, n
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      character(len=4096) :: msg
      character(len=256) :: linetext

      call find_unit(file_contaminant, 200)
      open(unit = file_contaminant, file="drutes.conf/ADE/contaminant.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
	print *, "missing drutes.conf/ADE/contaminant.conf file"
	ERROR STOP
      end if

      allocate(adepar(maxval(elements%material)))
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/water.conf/matrix.conf  &
	the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_contaminant, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
	errmsg=trim(msg))
      
      write(unit=msg, fmt=*) "HINT 1: Is the molecular diffusion value positive?", new_line("a"), &
			      "HINT 2 : Is the number of molecular diffusion values corresponding to the amount of layers?"
			      
      do i=1, ubound(adepar,1)
	call fileread(adepar(i)%difmol, file_contaminant, ranges=(/0.0_rkind, 1.0_rkind*huge(tmp)/), &
	  errmsg=trim(msg))
      end do
 
      write(unit=msg, fmt=*) "HINT 1: Are all values anisotropy defining anisotropical diffusivity positive? ", new_line("a"), &
	"HINT 2 : Have you defined enough values for anisotropy &
	(e.g. for 2D define angle and the maximal and minimal value of diffusivity, in total 3 values)?", new_line("a"),&
	"HINT 3: The number of lines with diffusivity has to correspond to the number of materials & 
	defined by your mesh"
      
      
      allocate(tmp_array(drutes_config%dimen + 1))
      do i=1, ubound(adepar,1)
	allocate(adepar(i)%diff_loc(drutes_config%dimen))
	call fileread(r=tmp_array, fileid=file_contaminant, ranges=(/0.0_rkind, huge(tmp)/), errmsg=trim(msg))
	adepar(i)%anisoangle = tmp_array(1)
	adepar(i)%diff_loc = tmp_array(2:drutes_config%dimen + 1)
	allocate(adepar(i)%diff(drutes_config%dimen, drutes_config%dimen))
	call set_tensor(adepar(i)%diff_loc, (/adepar(i)%anisoangle/), adepar(i)%diff)
      end do
      
      
      do i=1, ubound(adepar,1)
       call comment(file_contaminant)
       read(unit=file_contaminant, fmt=*, iostat=i_err) adepar(i)%cinit, adepar(i)%icondtype
       if (i_err /= 0) then
	 write(unit=terminal, fmt=*) "The number of lines for initial concentration has to be equal to the number of materials"
	 backspace(file_contaminant)
	 call comment(file_contaminant)
	 read(unit=file_contaminant, fmt=*, iostat=i_err) linetext
	 write(unit=terminal, fmt=*) "the following inappropriate line was specified in your config file", trim(linetext)
	 call file_error(file_contaminant)   
       end if
       select case (adepar(i)%icondtype)
	 case("cr", "ca")
	   CONTINUE
	 case default
	   write(unit=terminal, fmt=*) " Your initial concentration can be only:  "
	   write(unit=terminal, fmt=*) "  ca - for absolute concentration"
	   write(unit=terminal, fmt=*) "  cr - for relative concentration"
	   call file_error(file_contaminant)
       end select
	   
       end do
       
       call fileread(with_richards, file_contaminant, &
        errmsg="Specify [y/n] to define whether you prefer to compute convection from the Richards equation or &
         specify the convection directly here")
       
       if (.not. with_richards .and. drutes_config%name=="ADE_wr") then
	 write(unit=terminal, fmt=*) "You have specified ADE_wr =(advection dispersion reaction equation, where &
	       convection is computed from the Richards equation, but you want to specify convection here."
	 write(unit=terminal, fmt=*) "Solution: specify model type ADEstd instead of ADE_wr."
	 ERROR stop
       else if (with_richards .and. drutes_config%name=="ADEstd") then
	write(unit=terminal, fmt=*) "You have specified ADEstd =(advection dispersion reaction equation, where &
	       convection is specified here, but you want the convection to be computed from the solution &
	       of the Richards equation."
	 write(unit=terminal, fmt=*) "Solution: specify model type ADE_wr instead of ADEstd."
	 ERROR stop
       end if
         
       
       if (.not. with_richards) then
	 if (allocated(tmp_array)) deallocate(tmp_array)
	 allocate(tmp_array(2))
	 do i=1, maxval(elements%material)
	   call fileread(tmp_array, file_contaminant, errmsg="convection has to be defined for each layer")
	   adepar(i)%convection = tmp_array(1)
	   adepar(i)%water_cont = tmp_array(2)
	 end do
      end if
       
       call fileread(n, file_contaminant, ranges=(/0_ikind, 1_ikind*huge(n)/), & 
	 errmsg="the number of orders of reactions has to be positive  or zero")
       
       
       write(unit=msg, fmt=*) "You have requested ", n," different orders of reactions.", new_line("a"), &
	   "For each different order specify its reaction rate, the reaction rates are specified in a line."
       do i=1, ubound(adepar,1)
	 allocate(adepar(i)%lambda(n))
	 call fileread(adepar(i)%lambda, file_contaminant, errmsg=trim(msg))
       end do
       
       write(msg, *) "The number of lines for kinetic/equilibrium sorption parameters has to be equal to the number of materials"
       do i=1, ubound(adepar,1)
	 call fileread(adepar(i)%sorption%kinetic, file_contaminant, errmsg=trim(msg))
      end do
      
      do i=1, ubound(adepar,1)
	 call fileread(adepar(i)%sorption%name, file_contaminant, errmsg=trim(msg), options=(/"freund", "langmu"/))
      end do
      
      if (allocated(tmp_array)) deallocate(tmp_array)
      allocate(tmp_array(3))
      
      do i=1, ubound(adepar,1)
	call fileread(tmp_array, file_contaminant, errmsg=trim(msg), ranges=(/0.0_rkind, huge(0.0_rkind)/))
	adepar(i)%sorption%adsorb=tmp_array(1)
	adepar(i)%sorption%desorb=tmp_array(2)
	adepar(i)%sorption%third=tmp_array(3)
      end do
      
       
      call fileread(n, file_contaminant, ranges=(/1_ikind, huge(n)/), &
      errmsg="You have selected strange number of boundaries for ADE problem")
      
      if (with_richards) then
	i=2
      else
	i=1
      end if
      
      call readbcvals(unitW=file_contaminant, struct=pde(i)%bc, dimen=n, &
	dirname="drutes.conf/ADE/")
      
      
      

    end subroutine ADE_read				


  

end module ADE_reader