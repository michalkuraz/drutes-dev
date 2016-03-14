module ADE_reader
  public :: ADE_read
  
  contains

    subroutine ADE_read(pde_loc)
      use typy
      use globals
      use global_objs
      use core_tools
      use ADE_globals
      use readtools
      use pde_objs

      class(pde_str), intent(in out) :: pde_loc
      integer :: i_err
      integer(kind=ikind) :: i, n
      real(kind=rkind) :: tmp
      real(kind=rkind), dimension(:), allocatable :: tmp_array
      character(len=4096) :: msg
      character(len=256) :: linetext
      logical :: crelative = .false.
      logical :: with_richards_def
      
      

      pde_loc%problem_name(1) = "ADER_in_liquid"
      pde_loc%problem_name(2) = "Advection-dispersion-reaction equation (solute concentration in liquid phase)"

      pde_loc%solution_name(1) = "solute_concentration" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "c  [M/L^3]" !popisek grafu

      pde_loc%flux_name(1) = "conc_flux"  
      pde_loc%flux_name(2) = "concentration flux [M.L^{-2}.T^{-1}]"

      pde_loc%mass_name(1) = "conc_in_porous_medium"
      pde_loc%mass_name(2) = "concetration [M/L^3]"
      
      

      call find_unit(file_contaminant, 200)
      open(unit = file_contaminant, file="drutes.conf/ADE/contaminant.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
	print *, "missing drutes.conf/ADE/contaminant.conf file"
	ERROR STOP
      end if
     

      allocate(adepar(maxval(elements%material)))
      
      call fileread(n, file_contaminant)
      
      backspace(file_contaminant)
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/ADE/contaminant.conf  &
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
	   if (adepar(i)%icondtype == "cr") crelative=.true.
	   CONTINUE
	 case default
	   write(unit=terminal, fmt=*) " Your initial concentration can be only:  "
	   write(unit=terminal, fmt=*) "  ca - for absolute concentration"
	   write(unit=terminal, fmt=*) "  cr - for relative concentration"
	   call file_error(file_contaminant)
        end select
	   
       end do
       
       if (crelative) then
	 call  fileread(adepar(1)%cmax, file_contaminant, errmsg="Have you defined correct value for the maximal concentration?")
	 do i=2, ubound(adepar,1)
	  adepar(i)%cmax = adepar(1)%cmax
	 end do
       end if
       
  
       
       if (.not. crelative) then
	write(msg, *) "HINT 1: Specify [y/n] to define whether you prefer to compute convection from the Richards equation or", &
        " specify the convection directly here.", new_line("a"), &
         "   HINT 2: Since you don't use relative concentration for the initial condition, check whether you left the", &
        " line with cmax blank."
       else
	write(msg, *) "HINT 1: Specify [y/n] to define whether you prefer to compute convection from the Richards equation or &
         specify the convection directly here."
       end if
       
       call fileread(with_richards, file_contaminant, errmsg=trim(msg))
       
       select case(adjustl(trim(drutes_config%name)))
	 case("ADEstd", "ADEstd_kinsorb")
	   with_richards_def = .false.
	 case default
	   with_richards_def = .true.
	end select
       
       if (.not. with_richards .and. with_richards_def) then
	 write(unit=msg, fmt=*) "You have specified" , &
	      " convection to be  computed from the Richards equation, but you want to specify the convection here." ,&
	      new_line("a"), "   Solution: change model type setup in global.conf or comment out the convection values here."
	 call file_error(file_contaminant, msg)
       else if (with_richards .and. .not. with_richards_def) then
	write(unit=msg, fmt=*) "You have specified here to compute the convection fromt the Richards equation, ",  &
	       " but in your global config you have requested solving ADE equation only, and so convection has to be specified here.", &
	 new_line("a"),  "   Solution:  change model type setup in global.conf or supply the convection here."
	 call file_error(file_contaminant, msg)
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
      
      if (with_richards) then
	write(unit=msg, fmt=*) "HINT1: Have you commented out all lines with convection values? ", &
	"Since the convection is computed from the Richards equation.", new_line("a"), &
	"   HINT2: The number of orders of reactions has to be positive or zero."
      else
	write(unit=msg, fmt=*) "HINT1: Is your number of lines with convection definition corresponding with the number of layers?",&
	new_line("a"), &
	"   HINT2: The number of orders of reactions has to be positive or zero."
      end if
       
       call fileread(n, file_contaminant, ranges=(/0_ikind, 1_ikind*huge(n)/), & 
	 errmsg=trim(msg))
       
       
       write(unit=msg, fmt=*) "You have requested ", n," different orders of reactions.", new_line("a"), &
	   "For each different "
       
       
       write(unit=msg, fmt=*) "You have requested ", n," different orders of reactions.", new_line("a"), &
	   "Specify the list of orders of reactions,", new_line("a"), &
	   "!!!!EACH material requires its line!!!!",  new_line("a"), &
	   "    e.g. you want to use zero order and second order reaction for first material and ",  new_line("a"), & 
	   " first and second order reaction for second material, then specify the following lines", new_line("a"),new_line("a"),&
 	   " 0 2 ", new_line("a"), " 1 2 "
       do i=1, ubound(adepar,1)
	 allocate(adepar(i)%orders(n))
	 call fileread(adepar(i)%orders, file_contaminant, errmsg=trim(msg))
       end do
       
       write(unit=msg, fmt=*) "You have requested ", n," different orders of reactions.", new_line("a"), &
	   "Specify the list of rates of reactions,", new_line("a"), &
	   "!!!!EACH material requires its line!!!!",  new_line("a"), &
	   "    e.g. you want to use zero order with rate 0.2 and second order reaction with rate 0.6 for first material and ",&
	   new_line("a"), & 
	   "and you have specified the requested orders of reactions above",  new_line("a"), & 
	   " for the second material you would like to use only  opnly second order kinetics with rate 0.4 ",& 
	   new_line("a"), "then the following lines has to be specified",&
	   new_line("a"),new_line("a"),&
 	   " 0.2 0.6 ", new_line("a"), " 0.0 0.4 "
 	   
       do i=1, ubound(adepar,1)
	 allocate(adepar(i)%lambda(n))
	 call fileread(adepar(i)%lambda, file_contaminant, errmsg=trim(msg))
       end do   

       
       select case(adjustl(trim(drutes_config%name)))
	  case("ADEstd_kinsorb", "ADE_RE_std_kinsorb", "ADE_REstdH_kinsorb", "ADE_RE_rot_kinsorb", "ADE_RErotH_kinsorb")
	    adepar(:)%sorption%kinetic = .true.
	  case default
	    adepar(:)%sorption%kinetic = .false.
	end select
 
      
      write(msg, *) "HINT1: The number of lines for kinetic/equilibrium sorption parameters has to be", &
      " equal to the number of materials" , &
      new_line("a"), "   HINT2: The bulk density has to be positive value greater than zero."
      
      do i=1, ubound(adepar,1)
	 call fileread(adepar(i)%bd, file_contaminant, errmsg=trim(msg), &
	 ranges=(/epsilon(0.0_rkind), huge(0.0_rkind)/))
      end do
      
      
      write(msg, *) "HINT1: The number of lines for kinetic/equilibrium sorption parameters has to be equal &
	to the number of materials"&
       , new_line("a"), "   HINT2: Have you specified the sorption model name correctly?"
      
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
	if (.not. adepar(i)%sorption%kinetic .and. adepar(i)%sorption%desorb < 100*epsilon(1.0_rkind) .and. &
	abs(adepar(i)%sorption%adsorb) > 100*epsilon(1.0_rkind)) then
	  write(msg, *) "This is non-sence!! You have defined zero desorption rate but you want to use equilibrium sorption.", &
	  new_line("a"), &
	  "Divisions by zero are in general not accepted in a good society."
	  call file_error(file_contaminant, message=trim(msg))
	end if
      end do
      
      if (drutes_config%name == "ADEwrk" .or. drutes_config%name=="ADEstk") then
	do i=1, ubound(adepar,1)
	  call fileread(adepar(i)%csinit, file_contaminant, errmsg="Have you defined value for c_s_init")
	end do
      end if
      
      if (drutes_config%name == "ADEwrk" .or. drutes_config%name=="ADEstk") then
 	write(msg, *) "HINT 1: You have selected strange number of boundaries for ADE problem.", new_line("a"), &
		    "  HINT 2: Since you requested nonequilibrium sorption, then you should provide csinit for each material."
      else
 	write(msg, *) "HINT 1: You have selected strange number of boundaries for ADE problem.", new_line("a"), &
		    "   HINT 2: Since you requested equilibrium sorption, then comment or erase lines with csinit"
      end if
		    
      
      call fileread(n, file_contaminant, ranges=(/1_ikind, huge(n)/), &
      errmsg=trim(msg))
      
      
      call readbcvals(unitW=file_contaminant, struct=pde_loc%bc, dimen=n, &
	dirname="drutes.conf/ADE/")
      
      
      

    end subroutine ADE_read		
    
    subroutine ADEcs_read(pde_loc)
      use typy
      use pde_objs
          
      class(pde_str), intent(in out) :: pde_loc
      
      pde_loc%problem_name(1) = "ADER_in_solid"
      pde_loc%problem_name(2) = "Advection-dispersion-reaction equation (solute concentration adsorbed in solid phase)"

      pde_loc%solution_name(1) = "solute_concentration" !nazev vystupnich souboru
      pde_loc%solution_name(2) = "c  [M/M]" !popisek grafu

      pde_loc%flux_name(1) = "zero_flux"  
      pde_loc%flux_name(2) = "zero flux"

      pde_loc%mass_name(1) = "conc_in_solid_phase"
      pde_loc%mass_name(2) = "concetration [M/M]"
      
    
    end subroutine ADEcs_read


  

end module ADE_reader