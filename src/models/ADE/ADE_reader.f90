module ADE_reader
  public :: ADE_read, ADEcs_read
  
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
      integer(kind=ikind) :: i, n, equipos
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
      
      

      open(newunit = file_contaminant, file="drutes.conf/ADE/contaminant.conf", action="read", status="old", iostat=i_err)
      if (i_err /= 0) then
        print *, "missing drutes.conf/ADE/contaminant.conf file"
        ERROR STOP
      end if
     
      
      write(msg, fmt=*) "ERROR!! incorrect number of materials in drutes.conf/ADE/contaminant.conf  &
        the mesh defines", maxval(elements%material)  , "materials, and your input file defines", n, "material(s)."
	
     
      call fileread(n, file_contaminant, ranges=(/1_ikind*maxval(elements%material),1_ikind*maxval(elements%material)/),&
        errmsg=trim(msg))
        
      call read_sep(file_contaminant)
      
      write(unit=msg, fmt=*) "HINT 1: Is the molecular diffusion value positive?", new_line("a"), &
			      "HINT 2 : Is the number of molecular diffusion values corresponding to the amount of layers?"
			      
      do i=1, ubound(adepar,1)
        call fileread(adepar(i)%difmol, file_contaminant, ranges=(/0.0_rkind, 1.0_rkind*huge(tmp)/), &
          errmsg=trim(msg))
      end do
      
      
      
      call read_sep(file_contaminant)
 
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
      
      call read_sep(file_contaminant)
      
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
       
       call read_sep(file_contaminant)
       
       if (.not. crelative) then
          write(msg, *) "HINT 1: Specify [y/n] to define whether you prefer to compute convection from the Richards equation or", &
        " specify the convection directly here.", new_line("a"), &
         "   HINT 2: Since you don't use relative concentration for the initial condition, check whether you left the", &
        " line with cmax blank."
       else
          write(msg, *) "HINT 1: Specify [y/n] to define whether you prefer to compute convection from the Richards equation or &
         specify the convection directly here."
       end if
       
       
       
 
       
       write(unit=msg, fmt=*) " The number of orders of reactions has to be positive or zero."
         
       
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
         if (n > 0) then
          call fileread(adepar(i)%orders, file_contaminant, errmsg=trim(msg), checklen=.true.)
         end if
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
         if (n>0) then
           call fileread(adepar(i)%lambda, file_contaminant, errmsg=trim(msg), checklen=.true.)
         end if
       end do   

      call read_sep(file_contaminant)
      
      call fileread(n, file_contaminant, ranges=(/1_ikind, huge(n)/), &
      errmsg=trim(msg))
      
      
      call readbcvals(unitW=file_contaminant, struct=pde_loc%bc, dimen=n, &
        dirname="drutes.conf/ADE/")
      
      
      

    end subroutine ADE_read		
    
    subroutine ADEcs_read(pde_loc)
      use typy
      use globals
      use global_objs
      use core_tools
      use ADE_globals
      use readtools
      use pde_objs
          
      class(pde_str), intent(in out), dimension(:) :: pde_loc
      real(kind=rkind), dimension(:), allocatable :: tmp_array, tmp_array2
      integer(kind=ikind) :: i, j, D
      real(kind=rkind) :: tmp
      character(len=4096) :: msg
      character(len=4096) :: number
      integer :: filesorp, ierr
      character(len=3) :: breaker
   
      
      open(newunit=filesorp, file="drutes.conf/ADE/sorption.conf", action="read", status="old", iostat=ierr)
      
      if (ierr /= 0) then
        print *, "unable to open drutes.conf/ADE/sorption.conf, exiting......."
        ERROR STOP
      end if
      
      do i=2, ubound(pde_loc,1)
        write(number, fmt=*) i-1
        write(pde_loc(i)%problem_name(1), fmt=*) "ADER_in_solid_", cut(number)
        pde_loc(i)%problem_name(1) = cut(pde_loc(i)%problem_name(1))
        write(pde_loc(i)%problem_name(2), fmt=*)  "Advection-dispersion-reaction equation (solute concentration adsorbed in &
              solid phase) ", cut(number)
        pde_loc(i)%problem_name(2) = cut(pde_loc(i)%problem_name(2))

        write(pde_loc(i)%solution_name(1), fmt=*) "solute_concentration_", cut(number) !nazev vystupnich souboru
        pde_loc(i)%solution_name(1) = cut(pde_loc(i)%solution_name(1))
        pde_loc(i)%solution_name(2) = "c  [M/M]" !popisek grafu

        pde_loc(i)%flux_name(1) = "zero_flux"  
        pde_loc(i)%flux_name(2) = "zero flux"

        pde_loc(i)%mass_name(1) = "conc_in_solid_phase"
        pde_loc(i)%mass_name(2) = "concetration [M/M]"
      end do
      
      if (no_solids > 0) then
        allocate(sorption(ubound(adepar,1), no_solids))
      else
        allocate(sorption(ubound(adepar,1), 1))
      end if
      
      write(msg, *) "The number of lines for kinetic/equilibrium sorption parameters has to be", &
         " equal to the number of materials and the number of columns has to be equal to the number", & 
         " of solid media.", new_line("a"), &
         "   e.g. 2 layers, 3 solid media with different sorption properties", & 
         " (2 with kinetic and 1 with equilibrium at each layer the same):", new_line("a"), &
         "   y y n ", new_line("a"), &
         "   y y n ", new_line("a"), & 
         "   y y n "
       
      do i=1, ubound(adepar,1)
        call fileread(sorption(i,:)%kinetic, filesorp, errmsg=msg)
      end do
      
      call read_sep(filesorp)
      
      if (no_solids == 0) then
        do i=1, ubound(sorption,1)
          if (sorption(i,1)%kinetic) then
            write(msg, *)  "You have specified in drutes.conf/ADE/ADE.conf zero number of soild media", new_line("a"), &
            " to gether with [y] for use soprtion. It means you want to use equilibrium sorption only.", new_line("a"),  &
            " But further in your file drutes.conf/ADE/sorption.conf , you requsted kinetic sorption", new_line("a"), & 
            " You should know what you want in your life, I'm not a crystal ball." 
            call file_error(filesorp, msg)
          end if
        end do
      end if

      write(msg, fmt=*) "   Set ratio  of solid media, if single solid medium ratio=1,",  & 
          " if more solid media the sum of all ratios per line has to be equal 1.0.", new_line("a"), &
          " The ratio must be defined for each scattered solid medium, in your case", &
          " the number of the scattered solid medium counts:", ubound(sorption,2) 
      do i=1, ubound(adepar,1)
        call fileread(sorption(i,:)%ratio, filesorp, errmsg=msg, checklen=.true.)
        tmp = 0
        do j=1, ubound(sorption,2)
          tmp = tmp + sorption(i,j)%ratio
        end do 
        
        if (abs(tmp - 1.0) > 100*epsilon(tmp)) then
          call file_error(filesorp, msg)
        end if
      end do
      
      call read_sep(filesorp)
    
       msg = "Set bulk density (for each material) if your medium scattered into more media provide bulk densities in columns"
       do i=1, ubound(adepar,1)
         call fileread(sorption(i,:)%bd, filesorp, errmsg=msg, &
               ranges=(/epsilon(0.0_rkind), huge(0.0_rkind)/), checklen=.true.)
       end do
       
       call read_sep(filesorp)
        
      write(msg, fmt=*) "Set sorption model name", new_line("a"), &
      "-   langmu for Langmuir model", new_line("a"), &
      "-   freund for Freundlich model."
      do i=1, ubound(sorption,1)
       call fileread(sorption(i,:)%name, filesorp, errmsg=msg, options=(/"freund", "langmu"/))
      end do
  

      
      allocate(tmp_array(3))  
      
      write(msg, fmt=*) "If kinetic Freundlich type (both k_adsorb and k_desorb is supplied here as a positive value)", &
      " and n (for nonlinear Freundlich model, for linear set n=1). For Langmuir model set ka, kd and csmax.", &
      " If equilibrium Freundlich model specify just k_D and n value, for equilibrium Langmuir model specify ", &
      " k_D and csmax only.", &
      "If your solid medium is scattered into more media, then  provide the data per layer and ", &
       "per medium always on separate lines, layers are always separated by ---"
      do i=1, ubound(sorption,1)
        do j=1, ubound(sorption,2)
          if (sorption(i,j)%kinetic) then
            D = 3
          else
            D = 2
          end if
          call fileread(tmp_array(1:D), filesorp, errmsg=trim(msg), ranges=(/0.0_rkind, huge(0.0_rkind)/), checklen=.TRUE.)
          sorption(i,j)%adsorb=tmp_array(1)
          sorption(i,j)%desorb=tmp_array(2)
          if (sorption(i,j)%kinetic) then
            sorption(i,j)%desorb=tmp_array(2)
            sorption(i,j)%third=tmp_array(3)
          else
            sorption(i,j)%third=tmp_array(2)
          end if
        end do
        if (i<=ubound(sorption,1)) call fileread(breaker, filesorp, errmsg="Each layer data must be separated by --- ", &
        options=(/"---"/))
      end do
      

      write(msg, *) "Initial concentration at solid phase (supply value for each layer and each medium)", &
    " use lines for layers and columns for media. if equilibrium model used provide value 0 otherwise an error will be generated."
      do i=1, ubound(sorption,1)
        call comment(filesorp)
        call fileread(sorption(i,:)%csinit, filesorp, errmsg=msg, checklen=.true.)
        do j=1, ubound(sorption,2)
           if (.not. sorption(i,j)%kinetic .and. abs(sorption(i,j)%csinit) > epsilon(tmp)) then
             write(msg, *) "Incorrect definition of csinit. Layer:", i, "solid medium:", j, new_line("a"), &
               "   is equilibrium sorption, so set 0 for initial concentration on solid and define only initial", &
               " concentration for liquid in file drutes.conf/ADE/contaminant.conf"
               call file_error(filesorp, msg)
           end if
         end do
       end do
        
!
      
      
    
    end subroutine ADEcs_read
    



  

end module ADE_reader
